#! /bin/bash

#SBATCH --job-name=proc_ufs_derived    # Specify job name
#SBATCH --partition=hera    # Specify partition name
#SBATCH --ntasks=12             # Specify max. number of tasks to be invoked
#SBATCH --mem-per-cpu=8000     # Specify real memory required per CPU in MegaBytes
#SBATCH --time=06:30:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=fv3lam       # Charge resources on this project account
#SBATCH --output=out.proc_UFS_derived   # File name for standard output
#SBATCH --error=err.proc_UFS_derived     # File name for standard error output

# Original Author: Hannah Christensen
# Date: Jan 2021
# Purpose: pre-process ICON DYAMOND output to produce area-weighted coarse grained fields
# These fields will be further processed to produce SCM input files
# [obtain time averaged fields on 0.1x0.1deg lat lon grid]
# This script computes derived variables on the fine resolution grid
#    >> mixinf ratios
#    >> follow with proc_ICON_derived_2.ksh to calculate theta and theta_l
#

# Modified by Xia Sun
# Date: Sep 2023
# Purpose: pre-process UFS DYAMOND output to produce distance-weighted coarse grained fields
# These fields will be further processed to produce CCPP-SCM input files
# [obtain time averaged fields on 0.2x0.2deg lat lon grid]
# This script computes derived variables on the fine resolution grid
#    >> mixinf ratios
#    >> follow with proc_UFS_derived_2.ksh to calculate theta and theta_l

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
day_start=2016090615
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



## 1. we will compute weight files which can be used for all remappings - this is the time consuming part of the regridding.
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
if [[ ! -d $output_dir_nat ]]; then
  echo "create output_dir directory"
  mkdir ${output_dir_nat}
fi


#########################################
###>> mixing ratio computation
echo "mixing ratio computation"
current_time="${day_start}"
while [[ "${current_time}" -le "$day_end" ]]; do
  echo "${current_time}"

    # convert full resolution grb -> netcdf; calculate at high resolution, then coarse grain
    cdo -f nc4 copy ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.grib2 ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.nc
    cdo -f nc4 copy ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.grib2 ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.nc

    # loop over vertical levels
    for j in $(seq -f "%03g" ${lev_start} ${lev_end}) ; do
        echo "starting level ${j}"

        # extract vertical level
        ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.nc
        ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc
        ncks -d lev,${j} ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.nc

        # merge
        cdo merge  ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_rl_ri.${j}.nc

        # compute mixing ratio at full resolution
        cdo -expr,'rv=q/(1.0-q)' ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_rl_ri.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.nc
        cdo merge ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_rl_ri.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_all_r.${j}.nc 
        cdo -expr,'ql=clwmr/(1.0+rv+icmr+clwmr)' ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_all_r.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.nc
        cdo -expr,'qi=icmr/(1.0+rv+icmr+clwmr)' ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_all_r.${j}.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.nc

        # regrid
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.${resol_target}.tmp.nc
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.${resol_target}.tmp.nc
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.${resol_target}.tmp.nc
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.${resol_target}.tmp.nc
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.${resol_target}.tmp.nc
        cdo -f nc4 -P 4 -O remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.${resol_target}.tmp.nc

        # convert height to record dimension
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.${resol_target}.nc
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.${resol_target}.nc
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.${resol_target}.nc
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.${resol_target}.nc
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.${resol_target}.nc
        ncpdq -a lev,time ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.${resol_target}.nc
   
        # tidy a little
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.${resol_target}.tmp.nc
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.${resol_target}.tmp.nc
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.${resol_target}.tmp.nc
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.${resol_target}.tmp.nc
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.${resol_target}.tmp.nc
        rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.${resol_target}.tmp.nc      

        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.rv.${j}.nc
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ql.${j}.nc
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qi.${j}.nc
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_all_r.${j}.nc 
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.qv_rl_ri.${j}.nc

        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${j}.nc
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${j}.nc # needed for theta_l calculation
        rm ${output_dir_nat}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${j}.nc

    done

# merge all levels
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.tmp.nc
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.tmp.nc
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.tmp.nc
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.tmp.nc
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.tmp.nc
ncrcat ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.*.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.tmp.nc

# restore time as record dimension
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.nc
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.nc
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.tmp.nc  ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.nc
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.nc
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.tmp.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.nc
ncpdq -a time,lev ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.tmp.nc  ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.nc

# extract desired region
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.${region}.nc
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.${region}.nc
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.${region}.nc
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.${region}.nc
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.${region}.nc
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.${region}.nc


# And compute r_t for completeness (linearity holds)
cdo merge ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.nc
cdo -expr,'rt=rv+clwmr+icmr' ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rt.${resol_target}.nc
rm -r ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.nc 

cdo merge ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.rv.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.CLMR.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.ICMR.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.${region}.nc
cdo -expr,'rt=rv+clwmr+icmr' ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.rt.${resol_target}.${region}.nc
rm -r ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_r.${resol_target}.${region}.nc


# compute q_t
cdo merge ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.nc
cdo -expr,'qt=q+ql+qi' ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qt.${resol_target}.nc
rm -r ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.nc

cdo merge ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.SPFH.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.ql.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.qi.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.${region}.nc
cdo -expr,'qt=q+ql+qi' ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.qt.${resol_target}.${region}.nc
rm -r ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.all_q.${resol_target}.${region}.nc

# tidy: delete individual level files
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.rv.*.${resol_target}.nc
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ql.*.${resol_target}.nc
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.qi.*.${resol_target}.nc
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.SPFH.*.${resol_target}.nc
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.CLMR.*.${resol_target}.nc
rm ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.ICMR.*.${resol_target}.nc
current_time=$(date -d "${current_time:0:8} ${current_time:8} 3 hours" +%Y%m%d%H)
done
