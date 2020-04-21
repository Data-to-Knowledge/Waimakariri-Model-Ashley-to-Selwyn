# Waimakariri-Model-Ashley-to-Selwyn

The best access to information is the report in this REPO,

How to access the model and input data (Waimakariri Model Data Catalogue)
The Waimakariri Model is a set of stochastic realizations and was developed in python using a number of open source and custom built libraries. This model is not usefully managed within a graphical user interface (GUI) such as visual modflow as these GUIs are optimised for development and management of single models.  At present, someone who has significant experience with the Modflow suite, the flopy library, and Python generally, should be able to access and use the model. This will take a significant amount of time for that user as the scripts have not been optimised for external use. That said with future work it would be possible to develop a more polished and robust library for this model.  It is unlikely that this model will ever be accessible to those with at least moderate python skills.
Scripts
As previously stated, the bespoke library of scripts surrounding this model are an intrinsic component of the model.  It is worth noting that these scripts are a working library and as such the comments and documentation within the library were primarily created for the developers.  Some effort has been made on key scripts (see Scripts) to improve the documentation, but it has not reached the quality of a published library. 

The best access to this library is through the github repository: . As a back up to the github repository, a copy of all scripts have been included within the Waimakariri Model Data Catalog.  These scripts should be treated as a last resort and are not guaranteed to be up to date beyond 21-04-2020.

While a detailed description of these scripts are beyond this resource the following structure and description should help future users get started.  Key scripts and modules that are likely to be useful to future users have been highlighted and labeled as (key) below:

Analytical_solutions: some analytical solutions that were used in the model runs.
Data_catalog_support: some supporting scripts for the creation of this data catalogue.
env:
Sdp.py (key): this script allows the user to set key locations. These locations must be set prior to the model being used. See the script header for more details.
Env_paths.py: depreciated, but held for documentation purposes.
Flopy_mh: a modified version of the flopy. This was used as there was not time during the process to create a pull request to fix several bugs identified in flopy.
Model_Tools: a set of bespoke tools to support modelling, model specific implementation in extended_boundary_model_tools.py.
Other_functions: a set of misc functions that were part of the larger ECAN environment and which were needed for key modelling.
Waimak_extended_boundry:
combine_nsmc_results: historical scripts left for documentation purposes that were used to combine the data received from GNS into useful formats.
Extended_boundry_model_tools.py (key): This script contains the extended waimakariri model specific implementation of the ModelTools class. This class has support for many common actions needed when modelling on a regular grid such as translation from geographic to model coordinates and model visualisation.
model_and_NSMC_build: historical scripts left for documentation purposes that were used to setup and create the pre-optimisation model.
model_checks.py: historical scripts left for documentation purposes that were used to check the pre-optimisation model for build flaws.
model_runs: historical scripts left for documentation purposes that were used to run specific predictions to inform plan change 7.
Model_run_tools:
Data_extraction (key): these scripts exist to extract head, drn and sfr flow, and feature concentration data from the model. They also allow the compilation of multiple model runs into much easier to manage NetCDF4 files.  Much of the data extraction can be done both from raw model files and from the collated NetCDF4 files. 
Forward_quanity_support (key): These scripts act to support the suite of forward ground and surface water quantity predictions.  
Metadata_managment (key): these scripts manage model metadata (e.g. convergence) and broad areas within the model.
Model_bc_data (key): these scripts support access to all of the model boundary condition (BC) data and support the creation of new BC data to make model predictions.  This includes N load layers.
Model_setup (key): These scripts contain basic wrappers to simplify the creation of Modflow, MT3D, and Modpath models.  It also contains the scripting necessary for re-creating a model realisation from the parameter files.
modpath_source_zone_support (key): These scripts are key supports for the generation of the source zone delimitation work for PC7.
n_analysis_support (key): These scripts support the creation of transport simulations for PC7.
stream_depletion_support (key): These scripts support some of the stream depletion work that was started, but not completed for the PC7 process.  These scripts are rougher than the rest within the model_run_tools module.
non_model_work: historical scripts left for documentation purposes that were used to create non model datasets to support the modelling. These were mainly in support of the development of land surface recharge layers, end-member mixing analysis for validation and filtering purposes and expert judgment elicitation support.
Nsmc_exploration_results: historical scripts left for documentation purposes that were used to explore the model on receipt from GNS
python_enviroment_information (key): data for the development of the python environment.  See Python environment for more details.
support_for_writeup: historical scripts left for documentation purposes that were used to support the write up of technical reports in support of Plan Change 7.
Waimak_modeling_non_extended: historical scripts left for documentation purposes that were used for the model prior to 04/2017, e.g. before the southern boundary was extended from Christchruch to past the Selwyn River.  These scripts are very rough.
Python environment
The python environment is the set of libraries used to interpret and run scripts. In the past the python environment for the Waimakariri Model has been somewhat unstable. In order to combat this I have put together an Anaconda YML file, which lists all of the dependencies and should be able to load the environment for the foreseeable future. The biggest difficulty in maintaining/updating the environment is the current flopy instance. This has been copied into these scripts as flopy_mh as a number of edits were needed in order to fix bugs in the package and there was not time in our process to undergo a pull request.  Unfortunately this requires numpy 1.13.

In order to install the environment:

download miniconda
open the anaconda prompt and enter:
        conda env create -f {path to environment.yml}  
the yml referenced here is: ...\python_env\waimak_model_environment.yml

In addition I have made a copy of the python environment in ...\python_env\Miniconda2. This folder can be copied to the local drive in order to be useful. This is not the best option for creating and using the python environment but it is included here to provide additional redundancy in this key step.  This environment was made on a windows 10 machine and may or may not be supported for future windows versions.  

If you are to use this pre-prepared python environment:
move it to a high speed local drive (otherwise it is painfully slow)
the interpreter you want is ...\Miniconda2\envs\waimak_model\python.exe
