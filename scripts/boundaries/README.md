# Boundaries using PyBDY

This directory contains setup files to run pybdy to make NAARC boundaries.
These are the namelist_naarc.bdy which has pybdy settings, grid_name_map.json which pybdy uses to match variable names in the NetCDF files with variables it expects and the src_data_naarc.ncml which points pybdy at the parent data files which are aggregated during loading. These 3 files need to be edited to have paths and variables that match the desired inputs. They have been set up for use with GOSI10p3 data.

The additional files: gen_bdy.slurm and iter_gen_bdy_NAARC.sh are for running pybdy on JASMIN to allow for the high memory usage which is about 50Gb. 

Full documentation on using pybdy  including how to install can be found here: (https://noc-msm.github.io/pyBDY/) and (https://github.com/NOC-MSM/pyBDY?tab=readme-ov-file). 
