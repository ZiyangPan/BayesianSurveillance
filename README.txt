Supplementary code repository for the paper:
Multi-Signal Safety Surveillance with Bayesian Latent Factor Modeling and Bias Correction

Project overview
This repository contains supplementary R code for two components of the paper:
(1) a simulation study for sequential Bayesian surveillance under latent factor modeling and bias correction, and
(2) a real-data analysis workflow based on profile-level input.

The repository is organized to keep the code compact and reader-friendly. The simulation component can be run directly. The real-data component provides the full analysis workflow, but the restricted real data are not included in this repository.

If you use this code, please cite the corresponding paper.

Repository structure
project/
├── README.txt
├── simulation/
│   └── run_simulation_sequential.R
├── real-data/
│   └── run_realdata_sequential.R
└── helper/
    ├── utils_general.R
    ├── utils_simulation.R
    └── utils_real_data.R

File descriptions
README.txt
Provides an overview of the repository, input requirements, output files, and general usage notes.

simulation/run_simulation_sequential.R
Main script for the simulation study. This script generates one full Day-60 dataset and then performs sequential truncation analyses at days 5, 10, ..., 60.

real-data/run_realdata_sequential.R
Main script for the real-data analysis. This script loads a user-supplied .RData file containing a profile-level object for this analysis workflow, constructs fixed exposure-outcome pairs over the requested period range, and runs the two-stage Bayesian analysis across the requested periods.

helper/utils_general.R
General helper functions shared by the simulation and real-data workflows.

helper/utils_simulation.R
Simulation-specific helper functions, including data generation, Stage 1 MCMC, Stage 2 MCMC, and the one-truncation analysis wrapper.

helper/utils_real_data.R
Real-data-specific helper functions, including fixed pair selection, profile likelihood construction, Stage 1 MCMC, Stage 2 MCMC, and the one-period analysis wrapper.

Simulation
The simulation study script is:
simulation/run_simulation_sequential.R

The simulation script generates one full dataset over Day 1 to Day 60 and then performs sequential truncation analyses at D_trunc = 5, 10, ..., 60.

Real-data analysis
The real-data script is:
real-data/run_realdata_sequential.R

The real data used in this project are not distributed with this repository.

To run the real-data analysis, update:
input_file <- "..."
profile_object_name <- "..."

The .RData file should contain a profile-level object for this analysis workflow.

If you are interested in the real-data component, please contact the authors for questions regarding data access.

Inputs and outputs
The simulation script is self-contained and does not require external data files.

The real-data script requires a user-supplied .RData file containing a profile-level object.

The scripts save analysis outputs as .RData and .csv files in the specified output location.

Runtime
Simulation runtime depends on MCMC settings and hardware.
Real-data runtime depends on profile size, number of periods, and MCMC settings.

Default MCMC settings
The default MCMC settings provided in the scripts are example settings and may be adjusted as needed.

Contact and data access
The real data are not included in this repository.
Please contact the authors for questions regarding real-data access or related analysis details.
