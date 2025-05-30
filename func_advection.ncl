; func_advection.ncl

;===================================================================================================
; calculate horizontal advective forcing required for SCM
;
; from email from Filip Vana:
;
;   => forcing is instantaneous
;   => advection forcing is limited to advection along model levels
;
;===================================================================================================

function func_advection(in_data:float,u_data:float,v_data:float)

local earth_radius,pi,deg_to_rad,tmp_array,lat_array,lon_array,x_position,y_position,x_grad_horx_m,x_grad_hory_m,hor_tend_data

begin

    ; in_data        -> quantity to advect
    ; u_data         -> u winds on correct pressure levels
    ; v_data         -> v winds on correct pressure levels

    ; calculate horizontal advection from instantaneous fields
    ; print("calculate advection")

    ; a) calculate horizontal advection
    earth_radius = 6371.22*10^3      ; radius in m
    pi = 4*atan(1.0)
    deg_to_rad = 2.0*pi/360.0

    tmp_array = new(dimsizes(in_data),typeof(in_data))
    lat_array = conform(tmp_array,in_data&lat,2)
    lon_array = conform(tmp_array,in_data&lon,3)

    x_position = earth_radius*cos(deg_to_rad*lat_array)*deg_to_rad*lon_array
    y_position = earth_radius*deg_to_rad*lat_array

    x_grad_horx_m = center_finite_diff_n(in_data,x_position,False,0,3)
    x_grad_hory_m = center_finite_diff_n(in_data,y_position,False,0,2)

    hor_tend_data   = - (u_data*x_grad_horx_m+v_data*x_grad_hory_m)

    copy_VarCoords(in_data,hor_tend_data)

    delete([/earth_radius,pi,deg_to_rad,tmp_array,lat_array,lon_array/])
    delete([/x_position,y_position,x_grad_horx_m,x_grad_hory_m/])

  return(hor_tend_data)

end

function center_finite_diff_o6(in_data:float,x_data:float,mydim:integer,smoothme:logical)

local in_data,x_data,mydim,smoothme,dims_in,i_hor,in_data_tmp,x_data_tmp,out_data_tmp,in_data_rearr,out_data_rearr,x_data_rearr,mean_dx,in_data_sm

begin

    ; in_data        -> quantity to calc derivative of
    ; x_data         -> x distance
    ; smoothme       -> logical to smooth before calculating deriv

    in_data_tmp = in_data
    in_data_tmp!0   = "zero"
    in_data_tmp!1   = "one"
    in_data_tmp!2   = "two"
    in_data_tmp!3   = "three"

    x_data_tmp   = x_data
    out_data_tmp = new(dimsizes(in_data),typeof(in_data))
    copy_VarCoords(in_data_tmp,x_data_tmp)
    copy_VarCoords(in_data_tmp,out_data_tmp)


    ; assume 4D data
    if (mydim.eq.0) then
       in_data_rearr = in_data_tmp
       out_data_rearr = out_data_tmp
       x_data_rearr  =  x_data
    end if
    if (mydim.eq.1) then
       in_data_rearr  =  in_data_tmp(one|:,zero|:,two|:,three|:)
       out_data_rearr = out_data_tmp(one|:,zero|:,two|:,three|:)
       x_data_rearr   =   x_data_tmp(one|:,zero|:,two|:,three|:)
    end if
    if (mydim.eq.2) then
       in_data_rearr  =  in_data_tmp(two|:,zero|:,one|:,three|:)
       out_data_rearr = out_data_tmp(two|:,zero|:,one|:,three|:)
       x_data_rearr   =   x_data_tmp(two|:,zero|:,one|:,three|:)
    end if
    if (mydim.eq.3) then
       in_data_rearr  =  in_data_tmp(three|:,zero|:,one|:,two|:)
       out_data_rearr = out_data_tmp(three|:,zero|:,one|:,two|:)
       x_data_rearr   =   x_data_tmp(three|:,zero|:,one|:,two|:)
    end if

    dims_in = dimsizes(in_data_rearr)
    mean_dx = new(dims_in,typeof(in_data_rearr))

    ; smooth dim 0 if necessary
    if (smoothme) then
       in_data_sm = in_data_rearr
       do i_hor = 3,dims_in(0)-4 ;==============================================
          in_data_sm(i_hor,:,:,:) = (0.05*in_data_rearr(i_hor-3,:,:,:) + \
                                     0.30*in_data_rearr(i_hor-2,:,:,:) + \
                                     0.75*in_data_rearr(i_hor-1,:,:,:) + \
                                     1.00*in_data_rearr(i_hor  ,:,:,:) + \
                                     0.75*in_data_rearr(i_hor+1,:,:,:) + \
                                     0.30*in_data_rearr(i_hor+2,:,:,:) + \
                                     0.05*in_data_rearr(i_hor+3,:,:,:))/3.2
       end do
       i_hor = 2                  ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.5*in_data_rearr(i_hor-2,:,:,:) + \
                                  2.0*in_data_rearr(i_hor-1,:,:,:) + \
                                  3.0*in_data_rearr(i_hor  ,:,:,:) + \
                                  2.0*in_data_rearr(i_hor+1,:,:,:) + \
                                  0.5*in_data_rearr(i_hor+2,:,:,:))/8.0
       i_hor = dims_in(0)-3       ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.5*in_data_rearr(i_hor-2,:,:,:) + \
                                  2.0*in_data_rearr(i_hor-1,:,:,:) + \
                                  3.0*in_data_rearr(i_hor  ,:,:,:) + \
                                  2.0*in_data_rearr(i_hor+1,:,:,:) + \
                                  0.5*in_data_rearr(i_hor+2,:,:,:))/8.0
       i_hor = 1                  ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.25*in_data_rearr(i_hor-1,:,:,:) + \
                                  0.50*in_data_rearr(i_hor  ,:,:,:) + \
                                  0.25*in_data_rearr(i_hor+1,:,:,:))
       i_hor = dims_in(0)-2       ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.25*in_data_rearr(i_hor-1,:,:,:) + \
                                  0.50*in_data_rearr(i_hor  ,:,:,:) + \
                                  0.25*in_data_rearr(i_hor+1,:,:,:))
    else
       in_data_sm = in_data_rearr
    end if

    ;; calculate the derivative on smoothed (or not) data

    do i_hor = 3,dims_in(0)-4 ;==============================================
       out_data_rearr(i_hor,:,:,:) = (-(1.0/60.0)*in_data_sm(i_hor-3,:,:,:) + \
                                       (3.0/20.0)*in_data_sm(i_hor-2,:,:,:) - \
                                       ( 3.0/4.0)*in_data_sm(i_hor-1,:,:,:) + \
                                       (1.0/60.0)*in_data_sm(i_hor+3,:,:,:) - \
                                       (3.0/20.0)*in_data_sm(i_hor+2,:,:,:) +\
                                       ( 3.0/4.0)*in_data_sm(i_hor+1,:,:,:))
       mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
       out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
    end do

    i_hor = 2  ;=============================================================
    out_data_rearr(i_hor,:,:,:) =  (1.0/12.0)*in_data_sm(i_hor-2,:,:,:) \
                                  -( 2.0/3.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +( 2.0/3.0)*in_data_sm(i_hor+1,:,:,:) \
                                  -(1.0/12.0)*in_data_sm(i_hor+2,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    i_hor = dims_in(0)-3 ;===================================================
    out_data_rearr(i_hor,:,:,:) =  (1.0/12.0)*in_data_sm(i_hor-2,:,:,:) \
                                  -( 2.0/3.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +( 2.0/3.0)*in_data_sm(i_hor+1,:,:,:) \
                                  -(1.0/12.0)*in_data_sm(i_hor+2,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
 

    i_hor = 1   ;============================================================
    out_data_rearr(i_hor,:,:,:) = -(1.0/2.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +(1.0/2.0)*in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    i_hor = dims_in(0)-2   ;=================================================
    out_data_rearr(i_hor,:,:,:) = -(1.0/2.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +(1.0/2.0)*in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    i_hor = 0  ;=============================================================
    out_data_rearr(i_hor,:,:,:) = -in_data_sm(i_hor,:,:,:) \
                                  +in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor,:,:,:))
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    i_hor = dims_in(0)-1  ;==================================================
    out_data_rearr(i_hor,:,:,:) =  in_data_sm(i_hor,:,:,:) \
                                  -in_data_sm(i_hor-1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor,:,:,:)-x_data_rearr(i_hor-1,:,:,:))
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    out_data = out_data_rearr(zero|:,one|:,two|:,three|:)
    copy_VarCoords(in_data,out_data)

    return(out_data)

end

function center_finite_diff_o4(in_data:float,x_data:float,mydim:integer,smoothme:logical)

local in_data,x_data,mydim,smoothme,dims_in,i_hor,in_data_tmp,x_data_tmp,out_data_tmp,in_data_rearr,out_data_rearr,x_data_rearr,mean_dx,in_data_sm

begin

    ; in_data        -> quantity to calc derivative of
    ; x_data         -> x distance
    ; smoothme       -> logical to smooth before calculating deriv

    in_data_tmp = in_data
    in_data_tmp!0   = "zero"
    in_data_tmp!1   = "one"
    in_data_tmp!2   = "two"
    in_data_tmp!3   = "three"

    x_data_tmp   = x_data
    out_data_tmp = new(dimsizes(in_data),typeof(in_data))
    copy_VarCoords(in_data_tmp,x_data_tmp)
    copy_VarCoords(in_data_tmp,out_data_tmp)


    ; assume 4D data
    if (mydim.eq.0) then
       in_data_rearr = in_data_tmp
       out_data_rearr = out_data_tmp
       x_data_rearr  =  x_data
    end if
    if (mydim.eq.1) then
       in_data_rearr  =  in_data_tmp(one|:,zero|:,two|:,three|:)
       out_data_rearr = out_data_tmp(one|:,zero|:,two|:,three|:)
       x_data_rearr   =   x_data_tmp(one|:,zero|:,two|:,three|:)
    end if
    if (mydim.eq.2) then
       in_data_rearr  =  in_data_tmp(two|:,zero|:,one|:,three|:)
       out_data_rearr = out_data_tmp(two|:,zero|:,one|:,three|:)
       x_data_rearr   =   x_data_tmp(two|:,zero|:,one|:,three|:)
    end if
    if (mydim.eq.3) then
       in_data_rearr  =  in_data_tmp(three|:,zero|:,one|:,two|:)
       out_data_rearr = out_data_tmp(three|:,zero|:,one|:,two|:)
       x_data_rearr   =   x_data_tmp(three|:,zero|:,one|:,two|:)
    end if

    dims_in = dimsizes(in_data_rearr)
    mean_dx = new(dims_in,typeof(in_data_rearr))

    ; smooth dim 0 if necessary
    if (smoothme) then
       in_data_sm = in_data_rearr
       do i_hor = 2,dims_in(0)-3   ;==============================================
          in_data_sm(i_hor,:,:,:) = (0.5*in_data_rearr(i_hor-2,:,:,:) + \
                                     2.0*in_data_rearr(i_hor-1,:,:,:) + \
                                     3.0*in_data_rearr(i_hor  ,:,:,:) + \
                                     2.0*in_data_rearr(i_hor+1,:,:,:) + \
                                     0.5*in_data_rearr(i_hor+2,:,:,:))/8.0
       end do
       i_hor = 1                  ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.25*in_data_rearr(i_hor-1,:,:,:) + \
                                  0.50*in_data_rearr(i_hor  ,:,:,:) + \
                                  0.25*in_data_rearr(i_hor+1,:,:,:))
       i_hor = dims_in(0)-2       ;==============================================
       in_data_sm(i_hor,:,:,:) = (0.25*in_data_rearr(i_hor-1,:,:,:) + \
                                  0.50*in_data_rearr(i_hor  ,:,:,:) + \
                                  0.25*in_data_rearr(i_hor+1,:,:,:))
    else
       in_data_sm = in_data_rearr
    end if

    ;; calculate the derivative on smoothed (or not) data
    do i_hor = 2,dims_in(0)-3  ;=========================================
       out_data_rearr(i_hor,:,:,:) =  (1.0/12.0)*in_data_sm(i_hor-2,:,:,:) \
                                     -( 2.0/3.0)*in_data_sm(i_hor-1,:,:,:) \
                                     +( 2.0/3.0)*in_data_sm(i_hor+1,:,:,:) \
                                     -(1.0/12.0)*in_data_sm(i_hor+2,:,:,:)
       mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
       out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
    end do
 
    i_hor = 1   ;============================================================
    out_data_rearr(i_hor,:,:,:) = -(1.0/2.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +(1.0/2.0)*in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
 
 
    i_hor = dims_in(0)-2   ;=================================================
    out_data_rearr(i_hor,:,:,:) = -(1.0/2.0)*in_data_sm(i_hor-1,:,:,:) \
                                  +(1.0/2.0)*in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor-1,:,:,:))/2
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
 
 
    i_hor = 0  ;=============================================================
    out_data_rearr(i_hor,:,:,:) = -in_data_sm(i_hor,:,:,:) \
                                  +in_data_sm(i_hor+1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor+1,:,:,:)-x_data_rearr(i_hor,:,:,:))
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)
 
 
    i_hor = dims_in(0)-1  ;==================================================
    out_data_rearr(i_hor,:,:,:) =  in_data_sm(i_hor,:,:,:) \
                                  -in_data_sm(i_hor-1,:,:,:)
    mean_dx(i_hor,:,:,:)  = (x_data_rearr(i_hor,:,:,:)-x_data_rearr(i_hor-1,:,:,:))
    out_data_rearr(i_hor,:,:,:) = out_data_rearr(i_hor,:,:,:)/mean_dx(i_hor,:,:,:)


    out_data = out_data_rearr(zero|:,one|:,two|:,three|:)
    copy_VarCoords(in_data,out_data)

    return(out_data)
 
end


function hor_derivative(in_data:float,smoothme:logical)

local smoothme,earth_radius,pi,deg_to_rad,tmp_array,lat_array,lon_array,x_position,y_position,x_grad_horx_m,x_grad_hory_m

begin

    ; in_data        -> quantity to calc derivatives for

    earth_radius = 6371.22*10^3      ; radius in m
    pi = 4*atan(1.0)
    deg_to_rad = 2.0*pi/360.0

    tmp_array = new(dimsizes(in_data),typeof(in_data))
    lat_array = conform(tmp_array,in_data&lat,2)
    lon_array = conform(tmp_array,in_data&lon,3)

    x_position = tofloat(earth_radius*cos(deg_to_rad*lat_array)*deg_to_rad*lon_array)
    y_position = tofloat(earth_radius*deg_to_rad*lat_array)

    ; if smoothme = tru, smooth each dimension in turn
    x_grad_horx_m = center_finite_diff_o4(in_data,x_position,3,smoothme)
    x_grad_hory_m = center_finite_diff_o4(in_data,y_position,2,smoothme)

    copy_VarCoords(in_data,x_grad_horx_m)
    copy_VarCoords(in_data,x_grad_hory_m)

    delete([/earth_radius,pi,deg_to_rad,tmp_array,lat_array,lon_array/])
    delete([/x_position,y_position/])

    return([/x_grad_horx_m,x_grad_hory_m/])

end
