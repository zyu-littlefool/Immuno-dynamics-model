%% CD8T/TEFF replay — all 15 seed groups
% This self-contained script replays CD8T/TEFF from antigen signals cached
% by the main model. It does NOT rerun the main model or change cached Ag.
% All CD8T/TEFF numerical code is included below; no external CD8 helper is
% required.
%
% Run scope:
%   15 seed groups x 3 release modes = 45 cached runs.
%   m = 1, 4, ..., 43  -> burst
%   m = 2, 5, ..., 44  -> sustained
%   m = 3, 6, ..., 45  -> burst-sustained
%
% Actual signal used in this analysis:
%   AntigenSignalRaw, corresponding to Ag_pre + Ag_pre_2.
%
% Before running:
%   1. Run run_full to generate the main-model caches.
%   2. Save this file as CD8_replay_full.m.
%   3. Run this script in MATLAB.
%
% Expected input layout:
%   inputRoot/
%     run01/CD8_input_cache.mat
%     run02/CD8_input_cache.mat
%     ...
%     run45/CD8_input_cache.mat
%
% The script validates the exact run indices before performing the replay.

%% 1. Paths and fixed analysis settings
repoRoot = fileparts(mfilename('fullpath'));

inputRoot = fullfile(repoRoot, ...
    'output', 'full_15groups');

outputRoot = fullfile(repoRoot, ...
    'output', 'full_15groups_CD8');

nSeedGroups = 15;
nReleaseModes = 3;
expectedRunIndices = 1:45;

signalSelection = "Ag_pre_plus_Ag_pre_2";

%% 2. CD8T/TEFF parameters used in the actual replay
P = struct();
P.decayBaseTEFF = 0.23;
P.decayTEFFCan = 0.47;

P.initialCD8T = 28;
P.alphaCD8 = 1.5;
P.nitaCD8T = 0.1;
P.decayCD8T = 0.23;

P.AgCanI = 5e4;

P.antigenK = 0.01;

%% 3. Load the cached runs and sort them by their original m indices
[cacheFiles, runIndices] = discoverRunCaches(inputRoot);
nRuns = numel(cacheFiles);

if nRuns ~= numel(expectedRunIndices) || ...
        ~isequal(runIndices(:).', expectedRunIndices)
    error('Expected run indices %s, but found %s.', ...
        mat2str(expectedRunIndices), mat2str(runIndices(:).'));
end

inputs = cell(nRuns, 1);

for k = 1:nRuns
    inputs{k} = readAgPreCache(cacheFiles{k});
end

fprintf('Loaded %d antigen-signal caches for CD8T/TEFF replay.\n', nRuns);

%% 4. Create a new replay session without overwriting earlier results
if isempty(outputRoot)
    outputRoot = fullfile(inputRoot, 'CD8_replays');
end
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

[~, sessionTag] = fileparts(tempname(outputRoot));
sessionDir = fullfile(outputRoot, ['CD8_run_replay_', sessionTag]);
mkdir(sessionDir);

ReplayInfo = struct( ...
    'InputRoot', inputRoot, ...
    'SourceFiles', {cacheFiles}, ...
    'Parameters', P, ...
    'SignalSelection', signalSelection, ...
    'SeedGroups', nSeedGroups, ...
    'ReleaseModes', nReleaseModes, ...
    'CreatedAt', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
save(fullfile(sessionDir, 'replay_settings.mat'), 'ReplayInfo');
% Archive the driver itself when it is running as a saved .m script.
scriptPath = [mfilename('fullpath'), '.m'];
if isfile(scriptPath)
    copyfile(scriptPath, ...
        fullfile(sessionDir, 'CD8_replay_full_snapshot.txt'));
end

%% 5. Replay CD8T/TEFF independently for every cached Ag trajectory
rows = cell(nRuns, 1);

for k = 1:nRuns
    runIndex = runIndices(k);
    seedGroup = floor((runIndex - 1) / nReleaseModes) + 1;
    C = inputs{k};

    groupDir = fullfile(sessionDir, ...
        sprintf('run%02d_group%02d_%s', ...
        runIndex, seedGroup, char(C.mode)));

    % The numerical CD8T/TEFF calculation is in the local function below.
    [~, row] = CD8TEFF_from_main( ...
        C.Ag, C.AntigenSignalRaw, C.dt, C.sampleIndex, C.mode, ...
        groupDir, P, signalSelection);

    % Add provenance to the summary table only. The per-run MAT file keeps
    % exactly the variables listed in the save call below.
    row.RunIndex = runIndex;
    row.SeedGroup = seedGroup;
    row.ReleaseMode = C.mode;
    row.SourceCacheFile = string(cacheFiles{k});
    row.OutputDirectory = string(groupDir);
    rows{k} = row;

    fprintf('Completed m%02d | seed group %02d | %s\n', ...
        runIndex, seedGroup, C.mode);
end

%% 6. Export the run summary
Summary = vertcat(rows{:});
Summary = movevars(Summary, ...
    {'RunIndex', 'SeedGroup', 'ReleaseMode', ...
     'SourceCacheFile', 'OutputDirectory'}, 'Before', 1);
Summary = sortrows(Summary, 'RunIndex');

writetable(Summary, ...
    fullfile(sessionDir, ...
    sprintf('CD8T_TEFF_run_%drun_summary.csv', nRuns)));

disp(Summary);
fprintf(['CD8T/TEFF replay completed for all %d runs.\n' ...
    'Output directory: %s\n'], nRuns, sessionDir);

%% Local functions
function [files, runIndices] = discoverRunCaches(inputRoot)

entries = dir(fullfile(inputRoot, '*', 'CD8_input_cache.mat'));
files = cell(numel(entries), 1);
runIndices = NaN(numel(entries), 1);

for k = 1:numel(entries)
    files{k} = fullfile(entries(k).folder, entries(k).name);
    [~, runFolder] = fileparts(entries(k).folder);
    token = regexp(runFolder, 'run(\d+)(?:_|$)', ...
        'tokens', 'once', 'ignorecase');
    runIndices(k) = str2double(token{1});
end

[runIndices, order] = sort(runIndices);
files = files(order);
end

function C = readAgPreCache(filePath)
%READAGPRECACHE Load one full-resolution antigen-signal cache.

C = load(filePath);
C.mode = string(C.mode);
end

function [CD8Columns, Summary] = CD8TEFF_from_main( ...
    Ag, AntigenSignal, dt, sampleIndex, mode, outputDir, overrides, signalSource)
% Local CD8T/TEFF numerical implementation for the replay driver.
% Argument 8 labels the selected antigen signal: Ag_pre + Ag_pre_2.
% The single parameter structure P is supplied by the driver above. The full
% RAW antigen input is cached before numerical cleanup or CD8 integration.

if nargin < 7 || isempty(overrides)
    overrides = struct();
end
if nargin < 8 || isempty(signalSource)
    signalSource = "Caller-provided signal";
end
signalSource = string(signalSource);
P = overrides;

Ag = Ag(:);
AntigenSignalRaw = AntigenSignal(:);
nSteps = numel(Ag)-1;
sampleIndex = sampleIndex(:);
mode = string(mode);
outputDir = char(outputDir);
if ~isfolder(outputDir)
    mkdir(outputDir);
end
Step = (0:nSteps)';
Time = Step * dt;
CacheInfo = struct('SchemaVersion',2, ...
    'AntigenSource',char(signalSource), ...
    'SourceFolder',outputDir, ...
    'CreatedAt',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
% A MAT file preserves full precision, full resolution and tiny imaginary parts.
% This is ONLY the selected CD8 input, not the main model's two-signal cache.
cacheFile = fullfile(outputDir,'CD8_input_cache.mat');
save(cacheFile,'AntigenSignalRaw','Ag','dt','sampleIndex','mode', ...
    'Step','Time','CacheInfo','-v7');

% Numerical cleanup affects the CD8 input copy only; Ag is kept unchanged.
agReal = real(AntigenSignalRaw);
agImag = imag(AntigenSignalRaw);
negativeMask = agReal < 0;
imagMask = agImag ~= 0;
AntigenSignal = agReal;
AntigenSignal(negativeMask) = 0;
if any(negativeMask | imagMask)
    warning('CD8TEFF:AgRoundoffCleaned', ...
        'CD8 input cleanup: %d tiny negatives and %d tiny imaginary parts.', ...
        nnz(negativeMask),nnz(imagMask));
end

CI = zeros(nSteps+1,1);
CD8T = zeros(nSteps+1,1);
TEFF = zeros(nSteps+1,1);
TotalDecayTEFF = zeros(nSteps+1,1);
CD8T(1) = P.initialCD8T;
TotalDecayTEFF(1) = P.decayBaseTEFF + P.decayTEFFCan;

for i = 1:nSteps
    % Entry i+1 is the selected stimulus from main-model step i.
    stimulus = AntigenSignal(i+1);
    CI(i+1) = min(stimulus*P.AgCanI,1.2e4);
    denom = CD8T(i) + CI(i+1);
    if denom == 0
        proliferation = 0; % zero cells and zero stimulus
    else
        proliferation = P.alphaCD8*CD8T(i)*CI(i+1)/denom;
    end
    CD8T(i+1) = CD8T(i) + dt*(proliferation-P.nitaCD8T*CD8T(i));
    CD8T(i+1) = CD8T(i+1)*(1-dt*P.decayCD8T);
    TotalDecayTEFF(i+1) = P.decayBaseTEFF + ...
        P.decayTEFFCan*(1-stimulus/(stimulus+P.antigenK));
    % Preserve the original update order: TEFF uses OLD CD8T(i).
    TEFF(i+1) = TEFF(i) + ...
        dt*(P.nitaCD8T*CD8T(i)-TotalDecayTEFF(i+1)*TEFF(i));
end

CD8Columns = [CD8T(sampleIndex),TEFF(sampleIndex), ...
    AntigenSignal(sampleIndex),CI(sampleIndex),TotalDecayTEFF(sampleIndex)];
[peakCD8T,idxCD8T] = max(CD8T);
[peakTEFF,idxTEFF] = max(TEFF);
Summary = table(mode,P.decayBaseTEFF,P.decayTEFFCan, ...
    peakCD8T,Time(idxCD8T),peakTEFF,Time(idxTEFF), ...
    mean(TEFF(sampleIndex)),TEFF(end), ...
    trapz(Time(sampleIndex),TEFF(sampleIndex)),trapz(Time,TEFF), ...
    'VariableNames',{'Mode','DecayBaseTEFF','DecayTEFFCan', ...
    'PeakCD8T','PeakCD8TTime','PeakTEFF','PeakTEFFTime', ...
    'MeanTEFF','FinalTEFF','AUCTEFF','AUCFullResolution'});
Summary.AntigenSource = signalSource;
Trajectories = table(repmat(mode,numel(sampleIndex),1),Time(sampleIndex), ...
    Ag(sampleIndex),AntigenSignal(sampleIndex),CI(sampleIndex), ...
    CD8T(sampleIndex),TEFF(sampleIndex),TotalDecayTEFF(sampleIndex), ...
    'VariableNames',{'Mode','Time','Ag','AntigenSignal','CI','CD8T','TEFF', ...
    'TotalDecayTEFF'});
Parameters = P;
Parameters.dt = dt;
Parameters.nSteps = nSteps;
Parameters.sampleIndex = sampleIndex;
Parameters.antigenSource = CacheInfo.AntigenSource;
writetable(Summary,fullfile(outputDir,'CD8T_TEFF_fixed_summary.csv'));
writetable(Trajectories,fullfile(outputDir,'CD8T_TEFF_fixed_trajectories.csv'));
save(fullfile(outputDir,'CD8T_TEFF_fixed_results.mat'), ...
    'Summary','Trajectories','Parameters', ...
    'mode','Time','Ag','CD8T','TEFF');
end
