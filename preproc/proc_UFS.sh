#!/bin/bash
#! /bin/bash
#SBATCH --job-name=proc_ufs    # Specify job name
#SBATCH --partition=hera    # Specify partition name
#SBATCH --ntasks=12             # Specify max. number of tasks to be invoked
#SBATCH --mem-per-cpu=8000     # Specify real memory required per CPU in MegaBytes
#SBATCH --time=01:15:00        # Set a limit on the total run time
#SBATCH --mail-type=FAIL       # Notify user by email in case of job failure
#SBATCH --account=fv3lam       # Charge resources on this project account
#SBATCH --output=out.proc_UFS    # File name for standard output
#SBATCH --error=err.proc_UFS    # File name for standard error output

# Original Author: Hannah Christensen
# Date: Jan 2021
# Purpose: pre-process ICON DYAMOND output to produce area-weighted coarse grained fields
# These fields will be further processed to produce SCM inut files
# [obtain time averaged fields on 0.1x0.1deg lat lon grid]

# Modified by Xia Sun
# Date: Sep 2023
# Purpose: pre-process UFS LAM 3km-runs output to produce distance-weighted coarse grained fields
# These fields will be further processed to produce CCPP-SCM inut files
# [obtain time averaged fields on 0.2x0.2deg lat lon grid]
# Date: Mar 2024
# conservative interpolation in place using CDO with SCRIP grid description.
# SCRIPT grid description generated using curvilinear_to_SCRIP from ncl 

# Run using bash
#
module purge
module load gnu
module load intel/2023.2.0
module load netcdf/4.7.0
module load wgrib2
set -x # echo on
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
preproc_dir=${work_dir}/preproc

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

# clean existing files
rm -r ${output_dir_nat}/*
rm -r ${output_dir}/*

echo "compute CG files for ${day_start}"
echo "================================="

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
  echo "UFS_${resol_target}_grid_${region}_wghts.nc exists"
else
  # create subsetted weight file
  echo "UFS_${resol_target}_grid_${region}_wghts.nc"
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

###time invariant fields
if [[  -f ${output_dir}/${region}/${fixed_file}_${resol_target}_VGRID_${region}.nc ]]; then
  # regridded height data already exists
  echo "regridded height data exists exists"
else
  # create regridded height data
  # we have to create height data for UFS LAM from scratch
  # extract all variables with HGT:surface,SFCR:surface,LAND:surface
  wgrib2  ${fixed_file}.grib2 -match ":(HGT):" -grib ${fixed_file}.HGT.grib2.tmp
  wgrib2 ${fixed_file}.HGT.grib2.tmp -for_n 1:64 -grib ${fixed_file}.HGT.grib2
  cdo -f nc4 copy ${fixed_file}.HGT.grib2 ${fixed_file}_HGT.nc
  wgrib2 ${fixed_file}.grib2 -match '^(1335|1386|1533):' -netcdf ${fixed_file}_vgrid.nc

# create regrided height data
echo "create regridded height data for ${resol_target} grid"
cdo -P 4 -O -f nc remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${fixed_file}_vgrid.nc ${fixed_file}_${resol_target}_vgrid.nc
cdo -P 4 -O -f nc remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${fixed_file}_HGT.nc ${fixed_file}_${resol_target}_HGT.nc
cdo -merge ${fixed_file}_${resol_target}_vgrid.nc ${fixed_file}_${resol_target}_HGT.nc ${output_dir}/${fixed_file}_${resol_target}_VGRID.nc
echo "    extract vgrid desired region"
cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${fixed_file}_${resol_target}_VGRID.nc ${output_dir}/${region}/${fixed_file}_${resol_target}_VGRID_${region}.nc
fi
#clean up
rm -r ${fixed_file}.HGT.grib2.tmp
rm -r ${fixed_file}.HGT.grib2
rm -r ${fixed_file}_vgrid.nc
rm -r ${fixed_file}_HGT.nc

############################################
##>> 3D fields
##pres t w u v qv tot_qc_dia tot_qi_dia

 for vari in pres t w u v qv tot_qc_dia tot_qi_dia


for vari in PRES TMP VVEL DZDT UGRD VGRD SPFH CLMR ICMR
do
current_time="${day_start}"
while [[ "${current_time}" -le "$day_end" ]]; do
  echo "${current_time}"
  # extract variable first
  echo "extracting ${i} for ${vari} from grib2 file"
  wgrib2 ${input_dir}/${input_pfx}.${current_time}.${input_sfx}.grib2 -match ":(${vari}):" -grib ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.grib2.tmp
  wgrib2 ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.grib2.tmp -for_n 1:64 -grib ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.grib2
  rm -r ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.grib2.tmp


  echo "regridding ${i} for ${vari}"
  cdo -P 4 -O -f nc remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.grib2  ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.${resol_target}.nc
  echo "extract desired region"
  cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${vari}.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${vari}.${resol_target}.${region}.nc
  current_time=$(date -d "${current_time:0:8} ${current_time:8} 3 hours" +%Y%m%d%H)
done
done

############################################
##>> 2D fields atm2
# ICON: hfls(surface latent heat flux, positive upwards)  hfss(surface sensible heat flux, DEPHY sign convention: postive upwards) ps (surface pressure) pr(precipitation) clt(cloud cover)
# UFS: LHTFL:surface,SHTFL:surface, PRES:surface, APCP Total Precipitation (kh m-2) ,TCDC
# ICON: ts(surface temp) hus(qv_s) tauu(u momentum flux, ustar) tauv (v momentum flux)
# UFS:TMP:surface,SPFH:surface,UFLX:surface,VFLX:surface
current_time="${day_start}"
while [[ "${current_time}" -le "$day_end" ]]; do
  echo "${current_time}"

  grib_ids=(1391 1390 1334 1367 1423 1336 1338 1388 1389)
  vars=("LHTFL" "SHTFL" "ps" "APCP" "TCDC" "ts" "qv_s" "UFLX" "VFLX")
  for index in "${!grib_ids[@]}";do
    grib_id="${grib_ids[$index]}"
    var="${vars[$index]}"
    echo $var
    echo "extracting ${i} for ${var} from grib2 file"
    wgrib2  ${input_dir}/${input_pfx}.${current_time}.${input_sfx}.grib2 -for ${grib_id}:${grib_id} -grib ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${var}.grib2

    echo "regridding ${i} for ${var}"
    cdo -P 4 -O -f nc remap,${resol_target}_grid_domain.nc,UFS_${resol_target}_grid_domain_wghts.nc ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${var}.grib2  ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.nc

    echo "extract desired region"   
    cdo -sellonlatbox,${minlon},${maxlon},${minlat},${maxlat} ${output_dir}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.${region}.nc
    if [ ${grib_id} -eq 1423 ];then
      echo "removing extra lev dimension for TCDC"
      cdo vertsum ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.${region}.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.${region}.tmp.nc
      mv ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.${region}.tmp.nc ${output_dir}/${region}/${input_pfx}.${current_time}.${input_sfx}.${var}.${resol_target}.${region}.nc
    fi
  done
  current_time=$(date -d "${current_time:0:8} ${current_time:8} 3 hours" +%Y%m%d%H)
done
