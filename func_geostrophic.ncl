; func_geostrophic.ncl

function func_geostrophic(fout:file)



local fout,tv,tv_top,tv_bot,grav,R_dry,rpi,rday,rsiyea,rsiday,\
      romega,f_cori_vec,f_cori,smooth_me,\
      logpres,gradp,gradp_x,gradp_y,press_gradx_out,press_grady_out,vg_out,ug_out,\
      dsizes,dsizes2,tmparray,lsm_tmp,lsm,lsm_perm,d_size,pgx_out_lsm,pgy_out_lsm,\
      is_cyclic,guess_type,nscan,eps,relc, opt, rho_sfc


begin

   ;==================================================================================
   ; "geostrophic winds" : calculate forcing term from pressure gradient
   ;
   ; NB this file has been simplified cf. Cascade file as we are only considering ocean points
   ;  => our partial derivatives along model levels are at constant height
   ;
   ;  Take care if we want to move to include land points. Because ICON is not a true sigma co-ordinate system,
   ;   we lose some of the nice properties of derivatives along model levels that are used in
   ;   ECMWF documentation

    print("procedure geostrophic winds")

    ;=====================
    ; load required variables in from file
    pres_data_out           = fout->pressure_t
    lsm                     = fout->lsm

    t_data_out              = fout->temp_t
    qv_data_out             = fout->qv_t
    ql_data_out             = fout->ql_t
    qi_data_out             = fout->qi_t

    lev_out  = t_data_out&lev
    time_out = t_data_out&time
    n_lev_out = dimsizes(lev_out)

    ;=====================
    ; do calculation

    ; set up constants - following IFS convention
    grav   = 9.80665
    R_dry  = 287.0
    rpi    = 4*atan(1.0);
    rday   = 86400;
    rsiyea = 365.25*rday*2.0*rpi/6.283076;
    rsiday=rday/(1.0+rday/rsiyea);
    romega=2.0*rpi/rsiday;

    f_cori_vec=2*romega*sin(pres_data_out&lat*rpi/180);
    f_cori=conform_dims(dimsizes(pres_data_out),f_cori_vec,2)
    delete([/f_cori_vec/])

    ;=====================
    ; virtual temperature
    ; p. 17, Rogers and Yau

    qt_data = qv_data_out+ql_data_out+qi_data_out
    ; deduce specific mass of dry air
    qd_data = 1 - qt_data

    ; convert to m.r.
    rv_data_out = qv_data_out/qd_data

    epsilon = 0.622
    tv_top = t_data_out*(1 + rv_data_out/epsilon)
    tv_bot = 1 + rv_data_out
    tv = tv_top/tv_bot
    copy_VarCoords(t_data_out,tv)

    ;----
    ; (1/rho) * partial d p by dx at constant eta
    logpres = log(pres_data_out)
    ; HMC remove smoothing
    copy_VarCoords(pres_data_out,logpres)
    copy_VarAtts(pres_data_out,logpres)
    smooth_me = False
    gradp = hor_derivative(logpres,smooth_me)
    gradp_x = gradp[0]
    gradp_y = gradp[1]
    delete([/gradp,logpres/])

    press_gradx_out = R_dry*gradp_x*tv
    press_grady_out = R_dry*gradp_y*tv

    ;-- save estimate of rho_sfc for ustar computation
    rho_sfc = pres_data_out(:,n_lev_out-1,:,:)/(R_dry*tv(:,n_lev_out-1,:,:))
    rho_sfc!0="time"
    rho_sfc&time=pres_data_out&time
    rho_sfc!1="lat"
    rho_sfc&lat=pres_data_out&lat
    rho_sfc!2="lon"
    rho_sfc&lon=pres_data_out&lon

    delete([/tv,tv_top,tv_bot/])
    delete([/R_dry,rpi,rday,rsiyea,rsiday,romega/])
    delete([/gradp_x,gradp_y/])


    vg_out =  press_gradx_out/f_cori
    ug_out = -press_grady_out/f_cori

    vg_out_flt = tofloat(vg_out)
    ug_out_flt = tofloat(ug_out)

    copy_VarCoords(pres_data_out,ug_out_flt)
    copy_VarCoords(pres_data_out,vg_out_flt)
    copy_VarAtts(pres_data_out,ug_out_flt)
    copy_VarAtts(pres_data_out,vg_out_flt)
    ug_out_flt@long_name = "Geostrophic zonal wind"
    ug_out_flt@units     = "m/s"
    vg_out_flt@long_name = "Geostrophic meridional wind"
    vg_out_flt@units     = "m/s"

    add_to_file(fout,ug_out_flt   , "ug")
    add_to_file(fout,vg_out_flt   , "vg")

    delete([/ug_out,vg_out,ug_out_flt,vg_out_flt,f_cori/])

    return(rho_sfc)

end
