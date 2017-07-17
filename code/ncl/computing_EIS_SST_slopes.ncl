; This ncl program computes EIS and SST slopes using the time seris of low cloud cover, EIS and SST in the 20th century with 20 CMIP5 models. 
;**********************************************
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
  f1 = addfile("./low_cloud_time_series_1900_1999_CMIP5.nc","r")
  cl_ts = f1->low_cl_regional
  delete(f1)
; EIS
  f1 = addfile("./EIS_time_series_1900_1999_CMIP5.nc","r")
  EIS_ts = f1->EIS_regional
  delete(f1)
; SST 
  f1 = addfile("./SST_time_series_1900_1999_CMIP5.nc","r")
  SST_ts = f1->SST_regional
  delete(f1)
; regression analysis
  slope_cl_2_var = new((/2,5,18/),double)
  do model_no = 0,18-1
    do iregion=0,4
         x1=dtrend(EIS_ts(:,iregion,model_no),False)
         x2=dtrend(SST_ts(:,iregion,model_no),False)
         y =dtrend(cl_ts(:,iregion,model_no),False)
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

  SST_slope_output = addfile("./SST_EIS_slope_CMIP5.nc","c")
  fAtt               = True            ; assign file attributes
  fAtt@title         = "SST/EIS slope estimates"
  fAtt@creation_date = systemfunc ("date")
  fAtt@region_names = "Peru,Namibia,Australia,California,Canary"
  fAtt@CMIP5_model_names = "BCC-CSM1.1,CCSM4,CSIRO-Mk3.6,CanESM2,FGOALS-s2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM,MIROC-ESM-CHEM,MIROC5,MPI-ESM-LR,MRI-CGCM3,NorESM1-M"
  fileattdef( SST_slope_output, fAtt )
  SST_slope_CMIP5@unit="%/K"
  SST_slope_CMIP5!0="models"
  SST_slope_CMIP5&models = ispan(1,18,1)
  EIS_slope_CMIP5@unit="%/K"
  EIS_slope_CMIP5!0="models"
  EIS_slope_CMIP5&models = ispan(1,18,1)
  SST_slope_output ->SST_slope_CMIP5 = SST_slope_CMIP5
  SST_slope_output ->EIS_slope_CMIP5 = EIS_slope_CMIP5

 end

