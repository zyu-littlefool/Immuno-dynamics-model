## Overview

This repository contains:

- The main stochastic immune response model;
- Three antigen release modes: Burst, Sustained, and Burst-sustained;
- Example and full reproduction entry points;
- CD8 T-cell replay analysis using cached antigen signals;
- Reference outputs for checking the reproducibility of the example.

The example runs the three release modes for seed group 15,
corresponding to the original run indices 34–36.

The full reproduction runs all 15 seed groups and three release modes,
corresponding to 45 simulations in total.

## Repository Structure

    .
    ├── Code/
    │   ├── run.m
    │   ├── run_example.m
    │   ├── run_full.m
    │   ├── CD8_replay_example.m
    │   ├── CD8_replay_full.m
    │   ├── banlance_fast.m
    │   ├── decay.m
    │   ├── GC.m
    │   └── get_Shannon.m
    │
    ├── example_data/
    │   ├── README.md
    │   ├── BCR tree.R
    │   └── Sample_Data
    │
    ├── README.md
    ├── LICENSE
    └── DATA_LICENSE.md

## Software Requirements

Tested with:

- MATLAB: R2024a
- Parallel Computing Toolbox

The scripts use the MATLAB `local` cluster profile.
The number of workers can be adjusted in the entry-point scripts.
No internet connection is required for running the simulations after
the repository has been downloaded.

## Quick Start: Example Reproduction

Open MATLAB and run:
    cd Code
    run_example
This command runs only seed group 12:

| Run index | Release mode |
|---|---|
| 34 | Burst |
| 35 | Sustained |
| 36 | Burst-sustained |

The default simulation length is `t_max = 8000`.

Expected output directory:

    Code/output/example_seed15/

Each run folder contains the main CSV output, gene-related outputs,
and a `CD8_input_cache.mat` file for subsequent CD8 replay.

After completing `run_example`, run:

    CD8_replay_example

Expected input:

    Code/output/example_seed15/

Expected output:

    Code/output/example_seed15_CD8/


## Reference Data and Validation

The `example_data` directory contains reference outputs, not input data
required by the model.

To validate a fresh example run, compare the newly generated files with
the reference files described in:

    example_data/README.md

Numerical comparisons should focus on the specified model variables.

The following environment-dependent fields should not be compared:

- Absolute file paths;
- Timestamps;
- MATLAB session identifiers;
- Source cache paths;
- Other machine-specific metadata.

## Model Parameters

The default model parameters are defined in Supplementary Materials.

Additional parameters controlling GC, PC, and B-cell dynamics are defined
in the corresponding helper functions.

For reproduction of the reported results, use the default parameter values
without modification.

## License

The MATLAB source code in the `Code` directory is released under the
MIT License. See `LICENSE` for details.

The reference data in `example_data` are released under the Creative
Commons Attribution 4.0 International License (CC BY 4.0).
