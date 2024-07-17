# WESTPA/OpenMM Example Guide

## Overview

The `full_westpa_example` directory contains all the files necessary to run a Weighted Ensemble molecular dynamics simulation using WESTPA2 with the OpenMM MD engine.

Specifically, this example involves sampling the Abl1 kinase modeled with the Drude polarizable forcefield using a 3D progress coordinate.

The Drude topology for Abl1 was built using CHARMM-GUI's Drude Prepper.

## Just Propagation Example

Since there are many moving parts to a WESTPA simulation, I added a simplified example in the `just_propagation_example` directory. This contains all the files necessary to run a brute-force MD simulation of the Abl1 kinase domain with the Drude forcefield using OpenMM.

To run the simulation:

1. Change the "base_dir" path in the configuration file (`just_propagation_example/common_files/abl1_drude_config.json`).
2. Optionally, modify the simulation parameters/platform.
3. While in '<your_path>/gms_openmm_code_sample/just_propagation_example/common_files', run "prod_drude_no_westpa.py" (no arguments needed)

The purpose of `just_propagation_example` is to exemplify the propagation (MD) step of a WESTPA run. It does not include any analysis or resampling/splitting/etc.

## Full WESTPA Example

Conversely, `full_westpa_example` includes all components of a full WESTPA simulation. However, some adjustments are required to run it on a new system:

1. Update the path in `env.sh` to your conda environment containing WESTPA2, MDTraj, and OpenMM.
2. Update the `base_dir` path in `common_files/abl1_drude_config.json` to the root directory of your WESTPA2 simulation.
3. Adjust the SLURM submission file `slurm_zmq_8gpus.sh` to match your system's specifics, particularly `--n-workers=8` to match the number of GPUs you plan to use.
4. When submitting the SLURM job, do it from your simulation root folder

## Key Details

### Saving on Overhead

In WEMD, many parallel short simulations are run simultaneously. To save on overhead costs, it is prudent to load the OpenMM system from disk instead of rebuilding it at every segment.

This is the default setting in the configuration file - `prod_drude_westpa.py` loads `abl1_drude_system.xml`, built with the default parameters in the configuration file.

### Customizing OpenMM System Parameters

If you wish to change the OpenMM system parameters (such as adding/removing forces or changing the PME nonbonded cut-off), set `system_path` to null in the configuration file. This ensures that a new system will be built with your custom parameters.

### Saving Disk Space

To save disk space, the full WESTPA example only saves two frames per segment (`seg.pdb` and `parent.pdb`). Each consecutive segment is separated by a timescale of `prod_tau` (default is 100 ps).

Changing this requires changing the 'pcoord_len: 2' property in 'west.cfg' to the number of frames you want to save, and also updating 'xyz_save_freq' in 'abl1_drude_config.json'

### Why dt = 1fs?

We use an integration timestep of 1 fs (0.001 ps) since this simulation uses the Drude polarizable forcefield.

### Time Units in Configuration File

All values relating to time (`prod_tau`, `integration_ts`, `xyz_save_freq`, `log_save_freq`) in the configuration files are in picoseconds (ps).
