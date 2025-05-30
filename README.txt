README.txt

Julia Simonson, 30 May 2025

Removed unused, commented out code


Xia Sun, 13 September 2023

NCL files to produce DEPHY-standard SCM input for the MUMIP project using Unified Forecast System (UFS) outputs


Hannah Christensen, 5 September 2022

NCL files to produce DEPHY-standard SCM input for the MUMIP project

====

Changes for v2.0
- bug removed in surface fields:
       These fields were saved in DYAMOND every 15 min not every 3 hr 
       so have a different time coordinate to the other variables.
       Correction made in preprocessing scripts and safety flag added 
       to check time coordinate in func_read_hires.ncl

- all variables saved as floats not double
       Changes made in func_advtend.ncl and func_grostrophic.ncl

- New derived variables: theta_l, mixing ratios
       Can either be computed accurately using raw ICON data before coarsening,
       or approximately within ncl scripts using low-resolution state variables
        - To do former, run proc_ICON_derived.ksh and proc_ICON_derived_2.ksh
       prior to running ncl scripts. Then in ncl scripts, set 
       flag_accurate_r = True and flag_accurate_thetal = True
        - To do latter, in ncl scripts, set 
       flag_accurate_r = False and flag_accurate_thetal = False

- New output: advective tendencies
       Advective tendencies for theta_l and mixing ratios computed and archived.

- Assorted changes to run on Levante as opposed to Mistral
       Changes in GRIB variable interpretation file on Levante, led to
       changes in variable names

- Memory issues improved in preprocessing scripts
       For derived variables, compute fields one model level at a time.


====

Version 1.0
- perform initial coarsening using cdo
- consider small region over Indian Ocean
- only first few timesteps
- resultant files not yet tested on SCM


