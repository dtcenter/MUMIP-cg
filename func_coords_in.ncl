;===================================================================================================
; read in input-file coordinate systems (high-res data)
; if process multiple time stamps, time_data_tmp(0) ---->time_data_tmp(timestep_idx(t_count))
;===============================================================

function func_coords_in(generic_path:string,timestep_path:integer,timestep_idx:integer,t_step_in:integer,time_ref:double,\
                            filename:string,variable:string,\
                            CG_resol:string,UFS_exp:string,region:string)

local generic_path,timestep_path,timestep_idx,time_ref,filename,variable,\
      t_count,data_path_matr,data_path,in_file,t_step_in,\
      lon_in,lat_in,lev_in,\
      n_time_in,time_in,t_count


begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   READ IN DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   t_count = 0

   data_path_matr = (/generic_path,"/ufs_",CG_resol,"/",region,"/",filename,tostring(timestep_path(t_count)),tostring_with_format(t_step_in*timestep_idx(t_count), "%0.2i"),".",UFS_exp,".",variable,".",CG_resol,".",region,".nc"/)

   data_path      = str_concat(data_path_matr)

;; Input variable dimensions: (time, latitude, longitude) ;
   in_file        = addfile(data_path,"r")

;; Read in coordinate variables from this file
   lon_in = in_file->lon
   lat_in = in_file->lat
   lev_in = in_file->lev
 
   lev_in!0 = "nlev"

;;  Loop over required times and read in time coordinate
   n_time_in = dimsizes(timestep_idx)
   time_in = new(n_time_in,double)

   do t_count = 0,n_time_in-1
      print(t_count)
 
      data_path_matr = (/generic_path,"/ufs_",CG_resol,"/",region,"/",filename,tostring(timestep_path(t_count)),tostring_with_format(t_step_in*timestep_idx(t_count), "%0.2i"),".",UFS_exp,".",variable,".",CG_resol,".",region,".nc"/)

      data_path      = str_concat(data_path_matr)

      print(data_path)
      ;; Input variable dimensions: (time, latitude, longitude) ;
      in_file        = addfile(data_path,"r")
      time_data_tmp  = in_file->time

      ; combines ref time attribute with time to give YYYYMMDDHH
      mytime = cd_calendar(time_data_tmp(0),3)
      ; compute day and month offset
      diff_mmdd = floor((mytime-time_ref)/100)
      nmon = floor(diff_mmdd/100)
      nday = diff_mmdd-nmon*100

      ;; convert in time from hours to seconds
      ; days
      time_in(t_count) = nmon*31+nday
      ; hours
      time_in(t_count) = time_in(t_count)*24 + round(time_data_tmp(0),0) ; need to change if process multi time steps, see notes above  Xia. (Mar, 2024)
      ; seconds
      time_in(t_count) = 60*60*time_in(t_count)

      delete([/data_path_matr,data_path/])

   end do

   time_in@long_name = "time"
   time_in@units = "seconds since 2016-08-01 00:00:00"
   time_in@calendar = "proleptic_gregorian"
   time_in!0 = "time"
   time_in&time = time_in

   return([/lon_in , lat_in , lev_in , time_in/])

end
