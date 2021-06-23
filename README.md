# Table of contents:
- [Installation and Environment Setup](#installation-and-Environment-Setup)
  * [Conda](#conda)
  * [Python](#python)
  * [OpenMM](#openmm)
  * [Other Required Python Packages](#other-required-python-packages)
- [Prepare Input Structure](#prepare-input-structure)
- [Edit input_file.py](#edit-input_filepy)
- # [Edit submit_Exchange-min_replica.sh](#edit-submit_Exchange-min_replicash)
- [Job Submission](#job-submission)
- [File Descriptions](#file-descriptions)
<!-- toc -->
-----
# Installation and Environment Setup
## Conda
In order to use this code to run ensembles of constant pH MD simulations, you'll first need to make sure all needed dependencies are installed. This is easily accomplished with Conda. [Click here](https://www.anaconda.com/products/individual) for Conda download and installation information.

## Python
The code presented here was built and tested with Python v3.9.1. Other later versions will likely also work, but this has not been tested.

## OpenMM
This code uses the molecular dynamics (MD) engine openMM to run simulations. [Click here](http://docs.openmm.org/latest/userguide/application.html#installing-openmm) for openMM installation information.

**Note:** this code has been built and tested against openMM 7.5.0. Other later versions will likely also work, but this has not been tested.

## Other Required Python Packages:
There are many packages that are used that are part of the standard Python library, and are thus not described here. This code makes specific use of the following:
- pandas: 1.2.4

- **Note:** this code has been built and tested against the version described above, but as with other packages we don't anticipate using newer versions of pandas will be problematic.

## Slurm
This code has been built against the [Slurm job scheduler](https://slurm.schedmd.com/). As such, it is required. We have specifically built this code against Slurm 17.02.8, but we expect that this code will run using newer versions without problems.

## Constant pH Code Installation:
So long as all dependencies are installed, installing this code is as simple as downloading a copy of this Git repository, and copying needed input files into it. This is described below.

-----
# Prepare Input Structure

1. Create a pair of pdb/psf files with all the titrateble residues protonated (considered here: Glu, Asp, His, Cys, Lys).

   For CHARMM FF:
   - Use GLUP and ASPP patches for GLU and ASP
   - Use HSP topology for all HIS.
   - CYS involved in disulphide bonds the procedure will be omitted.

   **IMPORTANT**: The PDB structure data should start from the SECOND line. Similar to conventional MD OpenMM simulations, FIRST line should specify unit cell information.

   **IMPORTANT**: Make sure that the ionic strength of the resulting system is sufficient (e.g., ≥0.1 M). The method is heavily dependent on the overall concentration of ions in the system due to the PME calculation procedure (like any other constant-pH method).

2. Equilbrate structure as preferred.
----
# Edit input_file.py
To set up the constant pH (replica-exchnage) calculation change the input_file.py file:
1. Edit job and Slurm settings (comments in file provide specific details).

2. Specify a full path to your pdb/psf files (DO NOT forget the quotation marks).
    - If you want to set different geometries for different replicas put the paths to correspondong pdb files to pdb_state_files

3.  Put the PDB **residue numbers** of CYS involved in disulphide bond

4. Fill in the OpenMM parameters for conventional MD runs

5. Select a pH range and a pH-step between the replicas

6. Specify the names and values for residues with user-defined pKas (experimentally obtained or else)
       Make sure that the order of the names and values list is the same
7. If restart option is set to 'ON', lambda-list and state files for each replica should exist

8. Edit constant pH (replica-exchange) settings:
    - Specify the number of preliminary minimization steps

    - Set the number of iterations of MD - replica exchnage
    - 
    - Set the number of MD cycles before replica exchange attempt -(MD_nsteps_replicas) 
    - If needed, run preliminary MD-replica exchange cycles (prep_replicas = True). Useful when pH-dependent conformational shifts are expected but structure files are not available for all the states observed preliminary MD-lambda exchnage cycles (prep_lambdas = True). Useful when all pH states are determined but pKa values are far from default values
----
# Edit submit_Exchange-min_replica.sh
Specifically, just change the line `conda activate openmm_7_5_0` so that it activates the conda environment that contains openMM and other requisite software (described above). Also add another other commands here (such as `module load cuda`) if needed to setup environment.

----
# Job Submission 
To submit job (on Slurm submission node) execute `./submit_job.sh`. Note: do NOT `sbatch submit_job`. Specifically, this is just a simple bash script (not a Slurm job script) that submits all jobs to Slurm. It's quite light-weight and should run/finish instantly.

-----
# File Descriptions

• assign_prot_prob.py

  Script that reads the input pdb-file and creates the protonation probabilities for all
  the selected titratable residues

• createLambdaList.py

  Script that creates the lambda_list files with the protonation probabilities of all the
  titratable residues at all pH values selected

• definitions.py

  Input file that defines the differencies between the partial charges of the titratables
  residues in protonated and non-protonated states and input pKa values.

• Exchange-min-replica.py

  Run script

• fep_functions.py

  File with definitions of FEP functions (see OpenMM Alchemical Transformations)

• input_file.py

  MD simulation imput - the only file that NEEDS to be changed

• run_replica.py

  The script that runs individual replica at a particular pH as a one slurm job

• pHrex.py

  The functional file that defines CpH systems

• setup_pH_system.py

 The functional file containing the rules for altering the nonbonded interactions (see OpenMM Alchemical transformations)

-----
