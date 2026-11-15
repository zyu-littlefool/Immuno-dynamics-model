%% Full reproduction: all 15 seed groups and three release modes
% Runs all 45 original simulations, m = 1, 2, ..., 45.
% Save this file as run_full.m in the same directory as run.m.

repoRoot = fileparts(mfilename('fullpath'));

config = struct();
config.helperDir = repoRoot;
config.targetRoot = fullfile(repoRoot, 'output', 'full_15groups');
config.clusterProfile = 'local';
config.seedGroupIDs = 1:15;

% Use up to 45 workers, and respect the local cluster limit.
% Reducing workers changes runtime, not the selected random streams.
localCluster = parcluster(config.clusterProfile);
config.nWorkers = min(45, localCluster.NumWorkers);

Results_full = run(config);

fprintf('Full simulation completed. Results: %s\n', config.targetRoot);
