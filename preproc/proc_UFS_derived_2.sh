#!/bin/bash
#SBATCH --job-name=proc_ufs_derived2    # Specify job name
#SBATCH --partition=hera    # Specify partition name
#SBATCH --ntasks=12             # Specify max. number of tasks to be invoked
#SBATCH --mem-per-cpu=8000     # Specify real memory required per CPU in MegaBytes
#SBATCH --time=03:59:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=fv3lam       # Charge resources on this project account
#SBATCH --output=out.proc_UFS_derived_2    # File name for standard output
#SBATCH --error=job.proc_UFS_derived_2    # File name for standard error output

# Original Author: Hannah Christensen
# Date: Jan 2021
# Purpose: pre-process ICON DYAMOND output to produce area-weighted coarse grained fields
# These fields will be further processed to produce SCM input files
# [obtain time averaged fields on 0.1x0.1deg lat lon grid]
#
# This script computes derived variables on the fine resolution grid
#   >>  theta, theta_l
#   >> assumes proc_ICON_derived.ksh has already been run
#

# Modified by Xia Sun
# Date: Sep 2023
# Purpose: pre-process UFS DYAMOND output to produce distance-weighted coarse grained fields
# These fields will be further processed to produce CCPP-SCM input files
# [obtain time averaged fields on 0.2x0.2deg lat lon grid]
#
# This script computes derived variables on the fine resolution grid
#   >>  theta, theta_l
#   >> assumes proc_UFS_derived.sh has already been run

# Run using bash
#

module purge
module load gnu
module load intel/2023.2.0
module load netcdf/4.7.0
module load wgrib2
module load cdo
module load nco

## 0. User specified variables

resol_target=0.2
# Resolution of UFS LAM runs
ufs_resol="3km"
# Experiment name for UFS LAM runs
ufs_exp="DYAMOND_3km"
# UFS LAM outputs needs to be staaged in input_dir
work_dir="/scratch2/BMC/fv3lam/MUMIP/expt_dirs/cg-ufs"
input_dir=${work_dir}/ufs_${ufs_resol}
output_dir=${work_dir}/CG/ufs_${resol_target}
output_dir_nat=${work_dir}/CG/ufs_native
vgrid_dir=${work_dir}/ufs_input/${ufs_resol}
preproc_dir=${work_dir}/preproc_ufs

#file staged in preproc file to generate fixed files, such as interpolation weights, HGT fileds
fixed_file="srw.t12z.natlev.f006.mumip_io_3km"
input_pfx="srw.t12z.natlev"
input_sfx="mumip_io_3km"

# Required dates: NB. do not span month break.
day_start=2016090100
day_end=2016091021
# for looping over levels
lev_start=0
lev_end=63

domain_minlon=48.00
domain_maxlon=97.76
domain_minlat=-37.13
domain_maxlat=5.94

# Required region
region="IO"
minlon=51.0
maxlon=95.0
minlat=-35.0
maxlat=5.0


## 1. we will compute files which can be used for all remappings - this is the time consuming part of the regridding.
# first weight file (e.g., 3km_grid.nc)
if [[  -f ${resol_target}_grid.nc ]]; then
  # target grid file already exists
  echo "${resol_target}_grid.nc exists"
else
  # create target grid file
  echo "create ${resol_target}_grid.nc"
  cdo -O -f nc -topo,global_${resol_target} ${resol_target}_grid.nc
  cdo -sellonlatbox,${domain_minlon},${domain_maxlon},${domain_minlat},${domain_maxlat} ${resol_target}_grid.nc ${resol_target}_grid_domain.nc
fi

# second weight file (e.g., UFS_3km_grid_domain_wghts.nc), requires src_SCRIP.nc to generate
if [[  -f UFS_${resol_target}_grid_domain_wghts.nc ]]; then
  # weight file already exists
  echo "UFS_${resol_target}_grid_domain_wghts.nc exists"
else
  # create weight file
  # We don't have grid area or cell corner lats and lons for UFS LAM at this time, we are using Distance-weighted average remapping
  echo "create UFS_${resol_target}_grid_domain_wghts.nc"
  wgrib2  ${fixed_file}.grib2 -match '^(1335):' -netcdf ${fixed_file}.gh.nc
  cdo -P 16 --cellsearchmethod spherepart gencon,${resol_target}_grid_domain.nc -setgrid,src_SCRIP.nc ${fixed_file}.gh.nc UFS_${resol_target}_grid_domain_wghts.nc
fi


# repeat the above for a particular region (needed for theta computation)
# first weight file (e.g., 3km_grid_IO.nc)
if [[  -f ${resol_target}_grid_${region}.nc ]]; then
  echo "${resol_target}_grid_${region}.nc exists"
else
  # create subsetted target grid file
  echo "create ${resol_target}_grid_${region}.nc"
  cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${resol_target}_grid.nc ${resol_target}_grid_${region}.nc
fi

# second weight file (e.g., UFS_3km_grid_IO_wghts.nc), requires src_SCRIP.nc to generate
if [[  -f UFS_${resol_target}_grid_${region}_wghts.nc ]]; then
  echo "UFS_${resol_target}_grid_${region}_wghts.nc"
else
  # create subsetted weight file
  echo "create UFS_${resol_target}_grid_${region}_wghts.nc"
  # create grid description file UFS_grid
  cdo griddes ${resol_target}_grid_${region}.nc > UFS_grid
  cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${fixed_file}.gh.nc ${fixed_file}.gh.${region}.nc
  cdo -P 16 --cellsearchmethod spherepart gencon,${resol_target}_grid_${region}.nc -setgrid,src_SCRIP.${region}.nc ${fixed_file}.gh.${region}.nc UFS_${resol_target}_grid_${region}_wghts.nc
fi


## 2. Regrid using these weights.

# create output directories if needed
if [[ ! -d $output_dir ]]; then
  echo "create output directory"
  mkdir ${output_dir}
fi
if [[ ! -d $output_dir/$region ]]; then
  echo "create region output directory"
  mkdir ${output_dir}/$region
fi


##########################################
####>> theta computation
echo "theta computation"

current_time="${day_start}"

while [[ "${current_time}" -le "$day_end" ]]; do
  echo "${current_time}"

   # convert full resolution grb -> netcdf
   cdo -f nc4 copy ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.TMP.grib2 ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.nc
   cdo -f nc4 copy ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.PRES.grib2  ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES.nc 

    # loop over vertical levels
    for j in $(seq -f "%03g" ${lev_start} ${lev_end}) ; do

            # extract vertical level
            ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES.${j}.nc
            ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.${j}.nc

            # merge
            cdo merge ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES_TMP.${j}.nc

            # compute theta at full resolution
            cdo -expr,'theta=t*((100000/pres)^0.286)' ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES_TMP.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.nc

            # regrid
            cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.${resol_target}.tmp.nc

            # convert height to record dimension
            ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.${resol_target}.nc

            # tidy a little
            rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.${resol_target}.tmp.nc
            rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES.${j}.nc
            rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.PRES_TMP.${j}.nc 
    done

   # merge all levels
   ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${resol_target}.tmp.nc

   # restore time as record dimension
   ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${resol_target}.nc

   # extract desired region
   cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.theta.${resol_target}.${region}.nc

   # tidy: delete individual level files
   rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.theta.*.${resol_target}*

   current_time=$(date -d "${current_time:0:8} ${current_time:8} 3 hours" +%Y%m%d%H)
done

##########################################
####>> theta_l computation
echo "theta_l computation"
### already have native t, theta and rl ###

current_time="${day_start}"

while [[ "${current_time}" -le "$day_end" ]]; do
  echo "${current_time}"

   # convert full resolution grb -> netcdf
   cdo -f nc4 copy ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.grib2 ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.nc

    # loop over vertical levels
    for j in $(seq -f "%03g" ${lev_start} ${lev_end}) ; do

       # extract vertical level
       ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc

       # merge
       cdo merge ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta_TMP_CLMR.${j}.nc
       
       # compute theta_l at full resolution
       cdo -expr,'thetal=theta-(theta/t)*((2.501*10^6)/1005.7)*clwmr' ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta_TMP_CLMR.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.nc
       
       # regrid
       cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.${resol_target}.tmp.nc
       
       # convert height to record dimension
       ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.${resol_target}.nc
            
       # tidy a little
       rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.${resol_target}.tmp.nc
       rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.TMP.${j}.nc
       rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta.${j}.nc
       rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.theta_TMP_CLMR.${j}.nc
       rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc
       rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.thetal.${j}.nc
   done

    # merge all levels
    ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${resol_target}.tmp.nc
 
    # restore time as record dimension
    ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${resol_target}.nc
    
    # extract desired region
    cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.thetal.${resol_target}.${region}.nc

    # tidy: delete individual level files and tmp files
    rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.thetal.*.${resol_target}.nc

    current_time=$(date -d "${current_time:0:8} ${current_time:8} 3 hours" +%Y%m%d%H)
done
