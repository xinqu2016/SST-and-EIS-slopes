; This ncl program computes EIS and SST slopes using the time seris of low cloud cover, EIS and SST in the 20th century. 
;*******
; Concepts illustrated:
;   - Drawing a simple contour plot
;   - Making an axis logarithmic in a contour plot
;   - Changing the labels and tickmarks on a contour plot
;   - Creatint a main title
;   - Attaching coordinate arrays to a variable
;
;***********************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;************************************************
;===========================
begin
  
; low cloud cover
  f1 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/cl_1900_2099_monthly_plevel_5_regions_18_westpart_new.nc","r")
  f2 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/cl_1900_2099_monthly_plevel_5_regions_18_eastpart_new.nc","r")
  cl_ts = (dim_avg_n(f1->low_cl_regional(:,0:99,:,:),0)+dim_avg_n(f2->low_cl_regional(:,0:99,:,:),0))/2.0
  print(cl_mean)
  delete(f1)
  delete(f2)

; EIS
  f1 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/LTS_tas_T850_EIS_1900_2099_18_monthly_westpart_new.nc","r")
  f2 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/LTS_tas_T850_EIS_1900_2099_18_monthly_eastpart_new.nc","r")
  EIS =(f1->EIS_regional(0:2400-1,:,:)+f2->EIS_regional(0:2400-1,:,:))/2.0
  EIS_t = new((/12,200,5,18/),double)
  do i=0,200-1
    do j=0,11
      EIS_t(j,i,:,:) = EIS(j+i*12,:,:)
    end do
  end do
  EIS_ts = dim_avg_n(EIS_t(:,0:99,:,:),0)
  delete(EIS)
  delete(EIS_t)
  delete(f1)
  delete(f2)

; SST 
  f1 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/SST_1900_2099_18_monthly_westpart_new_free_mode.nc","r")
  f2 = addfile("/Users/xinqu/Projects/gig2_backup/low_cloud_IPCC_AR5/SST_1900_2099_18_monthly_eastpart_new_free_mode.nc","r")
  SST = (f1->SST_regional(0:2400-1,:,:)+f2->SST_regional(0:2400-1,:,:))/2.0

  SST_t = new((/12,200,5,18/),double)
  do i=0,200-1
    do j=0,11
      SST_t(j,i,:,:)  = SST(j+i*12,:,:)
    end do
  end do
  SST_ts = dim_avg_n(SST_t(:,0:99,:,:),0)
  delete(f1)
  delete(f2)
  delete(SST)

; regression analysis
  slope_cl_2_var = new((/2,5,18/),double)
  do model_no = 0,18-1
    do iregion=0,4
         x1=dtrend(EIS_ts(1900-1900:1999-1900,iregion,model_no),False)
         x2=dtrend(SST_ts(1900-1900:1999-1900,iregion,model_no),False)
         y =dtrend(cl_ts(1900-1900:1999-1900,iregion,model_no),False)

         X = new((/3,100/),double) 
         X(0,:) = 1
         X(1,:) = x1
         X(2,:) = x2
         b = reg_multlin(y,X,False)
         slope_cl_2_var(:,iregion,model_no) = b(1:2)

         delete(x1)
         delete(x2)
         delete(X)
         delete(y)
    end do
  end do
  EIS_slope_CMIP5 = dim_avg_n(slope_cl_2_var(0,:,:),0)
  SST_slope_CMIP5 = dim_avg_n(slope_cl_2_var(1,:,:),0)

 end
