##########################################################################################
#                                                                                        #
#                                                                                        #
#                                  CALCULATION SETUP                                     #
#                                                                                        #
#                                                                                        #
##########################################################################################

1. Create a pair of pdb/psf files with all the titrateble residues protonated.

   IMPORTANT: The PDB structre data should start from the SECOND line

   For CHARMM FF:
   Use GLUP and ASPP patiches for GLU and ASP; set HSP name for all HIS
   For CYS involved in disulphide bonds the procedure will be omitted

   IMPORTANT: Make sure that a ionic strength of the resulting system in sufficient (≥0.1 M).
              The method is heavily dependent on the overall concentration of ions in the system 
              due to PME calculation procedure (like any other constant-pH method) 

2. Perform standard (and prefered) equilibration procedure

3. To set up the constant pH replica/lambda-exchnage calculation change the input_file.py file:
    - a) Specify a full path to your pdb/psf files (DO NOT forget the quotation marks)
       If you want to set different geometries for different replicas - put the paths to correspondong pdb files to pdb_state_files
    - b) Put the PDB numbers of CYS residues involved in disulphide bond (Temporary)
    - c) Fill in the OpenMM parameters for conventional MD runs
    - d) Set the soft-core potential parameters for alchemical scaling of nonbonded forces for titratable  protons
    - e) Select a pH range and a pH-step between the replicas
    - f) Specify the names and values for residues with user-defined pKas (experimentally obtained or else)
       Make sure that the order of the names and values list is the same
    - g) If restart option is set to 'ON', lambda-list files should exist
    - h) The constant pJ replica/lambda exchnage procedure:
       Specify the number of preliminary minimization steps
       Set the percent amount of randomly titratable chosen residues for lambda-exchnage attempt (default is 0.1 - 10%)
       Set the number of iterations of MD - replica exchnage - MD lambda exchnage cycle
       Set the number of MD steps before replica exchange attempt (MD_nsteps_replicas) and lambda-exchnage attempt (MD_nsteps_lambdas)

       If needed run
          preliminary MD-replica exchange cycles (prep_replicas = True). Useful when pH-dependent conformational shifts are expected but structure files are not available for all the states observed
          preliminary MD-lambda exchnage cycles (prep_lambdas = True). Useful when all pH states are determined but pKa values are far from default values

4. To set slurm job run "sbatch submit_job.sh"
       

