%% Reproducibility example: seed group 12, all three release modes
% Runs the original simulations m = 34, 35, and 36 only.
% Save this file as run_example.m in the same directory as run.m.

repoRoot = fileparts(mfilename('fullpath'));

config = struct();
config.helperDir = repoRoot;
config.targetRoot = fullfile(repoRoot, 'output', 'example_seed15');
config.clusterProfile = 'local';
config.seedGroupIDs = 12;

% Use at most three workers, and respect the local cluster limit.
localCluster = parcluster(config.clusterProfile);
config.nWorkers = min(3, localCluster.NumWorkers);

Results_example = run(config);

fprintf('Example completed. Results: %s\n', config.targetRoot);
fprintf('Expected run folders: run34, run35, and run36.\n');
