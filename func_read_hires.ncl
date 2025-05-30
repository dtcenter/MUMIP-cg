;===================================================================================================
; read in atmospheric data from high-resolution files
;
; performs temporal interpolation if flag_interp = .TRUE.
; reads in only level indicated by lev_count
; if process multiple time stamps, in_file->$variable_name$(0,lev_count,:,:) ---->in_file->$variable_name$(timestep_idx(0),lev_count,:,:)
;===============================================================

function func_read_hires(generic_path:string,timestep_path:integer,timestep_idx:integer,t_step_in:integer,time_in:double,\
                            filename:string,variable_file:string,variable_name:string,\
                            CG_resol:string,UFS_exp:string,region:string,lev_count:integer)


local generic_path,timestep_path,timestep_idx,t_step_in,time_in,filename,variable_file,variable_name,\
      CG_resol,UFS_exp,region,lev_count,no_time_slices,t_count,data_path_matr,data_path,in_file,\
      in_data_tmp,var_dims,data_size,data_size_new,in_data,i_dim


begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   READ IN DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; 1. work out how many time slices we have in total
  no_time_slices = dimsizes(timestep_idx)

  print(no_time_slices)

  do t_count = 0,dimsizes(timestep_path)-1

   data_path_matr = (/generic_path,"/ufs_",CG_resol,"/",region,"/",filename,tostring(timestep_path(t_count)),tostring_with_format(t_step_in*timestep_idx(t_count), "%0.2i"),".",UFS_exp,".",variable_file,".",CG_resol,".",region,".nc"/)
     data_path      = str_concat(data_path_matr)
     print(data_path)

  ;  Input variable dimensions: (time, latitude, longitude) ;
     in_file        = addfile(data_path,"r")

     ; data check
     if (t_count.eq.0) then
        tmptmp = in_file->$variable_name$
        ; check dimension time
        tmpsize = dimsizes(tmptmp)
        ;print(tmpsize)
        if (tmpsize(0).ne.8)
           print(" ERROR: UNEXPECTED NUMBER OF INPUT TIME SLICES")
           tmp = "nope"
        end if
        delete([/tmptmp, tmpsize/])
     end if

     if (lev_count.eq.999) then
        ; 2D field: only one level available
        in_data_tmp    = in_file->$variable_name$(0,:,:)  ; need to change if process multi time steps, see notes above. Xia. (Mar, 2024)
     else
        ; 3D field
        in_data_tmp    = in_file->$variable_name$(0,lev_count,:,:) ; need to change if process multi time steps, see notes above. Xia. (Mar, 2024)
     end if

     ; set up the new empty array into which we will put our multiple time data
     if (t_count.eq.0) then

        var_dims = getvardims(in_data_tmp)
        data_size = dimsizes(in_data_tmp)
        data_size_new = array_append_record(no_time_slices,data_size,0)
        delete([/data_size/])
        data_size = data_size_new
        delete([/data_size_new/])
        in_data = new(data_size,float)
        in_data!0    = "time"
        do i_dim = 0,dimsizes(var_dims)-1
           in_data!(i_dim+1)           = in_data_tmp!(i_dim)
           in_data&$in_data!(i_dim+1)$ = in_data_tmp&$in_data_tmp!(i_dim)$
        end do
        delete([/var_dims,i_dim,data_size/])

        copy_VarAtts(in_data_tmp,in_data)
     end if

     ; put time slice into in_data
     data_size = dimsizes(in_data)

     if (dimsizes(data_size).eq.3) then     ; time,      lat, lon
        in_data(t_count,:,:)   = (/in_data_tmp(lat|:,lon|:)/)
     else                                   ; time, lev, lat, lon
        in_data(t_count,:,:,:) = (/in_data_tmp(lev|:,lat|:,lon|:)/)
     end if
     
     delete([/data_path_matr,data_path,in_file,data_size,in_data_tmp/])
  
  end do

  in_data&time = time_in

  return(in_data)

end
