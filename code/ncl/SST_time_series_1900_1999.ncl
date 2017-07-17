; This ncl program computes annual-mean SST time series based on 20th-century simulations with 20 CMIP5 models
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

 SST_regional = new((/100, 5,18/),double) 

 diri_init  = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo")
 model_index = (/2,5,7,8,10,11,12,13,15,16,17,19,22,23,24,25,28,29/)
 diri = diri_init(model_index)

 start_year    = (/1850,1850,1850,1850,1850,1860,1861,1896,1850,1885,1885,1850,1850,1850,1850,1850,1850,1850/)
 length1 = 1999-1900+1

 do model_no = 0,17;
; historical
  file_ts_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/hist/ts*.nc")

  f_ts_hist_0    = addfile(file_ts_hist(0), "r")   ; note the "s" of addfile
  f_ts_hist      = addfiles(file_ts_hist, "r")   ; note the "s" of addfile
  lat  = f_ts_hist_0->lat
  lon  = f_ts_hist_0->lon
  SST_3d = new((/12,100,dimsizes(lat),dimsizes(lon)/),double)

  do i=0,length1-1
    print ("model_no="+model_no +"    year="+i)
    do imonth=0,11
        ts_djfmam_mean = f_ts_hist[:]->ts((1900-start_year(model_no)+i)*12+imonth,:,:);
        SST_3d(imonth,i,:,:) = ts_djfmam_mean 
        delete(ts_djfmam_mean)
   end do
  end do;
  delete(file_ts_hist)
  delete(f_ts_hist)
  delete(f_ts_hist_0)

  file_sftlf = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/sftlf*")
  f_cl_sftlf_0   = addfile(file_sftlf(0), "r")   ; note the "s" of addfile
  sftlf = f_cl_sftlf_0->sftlf
  if (diri(model_no).eq."MIROC-ESM" .or. diri(model_no).eq."MIROC-ESM-CHEM")
     sftlf  = sftlf*100
  end if

  lat_2d  = new((/dimsizes(lat),dimsizes(lon)/),double)
  lon_2d  = new((/dimsizes(lat),dimsizes(lon)/),double)
  do k=0,dimsizes(lon)-1
     lat_2d(:,k)=lat(:)
  end do
  do k=0,dimsizes(lat)-1
     lon_2d(k,:)=lon(:)
  end do

  do iyear=0,length1-1
;  Peru
     EIS0   = dim_avg_n(SST_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lon_2d.ge.360-110 .and. lon_2d.le.360-70 .and. lat_2d.ge.-30 .and. lat_2d.le.-10 .and. sftlf.lt.0.01, EIS0,-999)
     weight = where(lon_2d.ge.360-110 .and. lon_2d.le.360-70 .and. lat_2d.ge.-30 .and. lat_2d.le.-10 .and. sftlf.lt.0.01, weight, -999)
     weight@_FillValue = -999
     EIS0@_FillValue = -999
     SST_regional(iyear,0,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)

; Namibia
     EIS0   = dim_avg_n(SST_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .ge. 360-25. .or. lon_2d.le.15) .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .ge. 360-25. .or. lon_2d.le.15) .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     SST_regional(iyear,1,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)

; Australia
     EIS0   = dim_avg_n(SST_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     SST_regional(iyear,2,model_no) = sum(EIS0*weight)/sum(weight)  
     delete(EIS0)
     delete(weight); 

; California
     EIS0   = dim_avg_n(SST_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and.  sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     SST_regional(iyear,3,model_no) = sum(EIS0*weight)/sum(weight)  
     delete(EIS0)
     delete(weight)

; Canary
     EIS0   = dim_avg_n(SST_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     SST_regional(iyear,4,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)

  end do

  delete(sftlf)
  delete(lat)
  delete(lon)
  delete(lat_2d)
  delete(lon_2d)
  delete(SST_3d)

 end do;
 f_output = addfile("SST_time_series_1900_1999_CMIP5.nc","w");
 fAtt               = True            ; assign file attributes
 fAtt@title         = "annual-mean SST time series in 18 CMIP5 models during the period of 1900-1999"
 fAtt@creation_date = systemfunc ("date")
 fAtt@region_names = "Peru,Namibia,Australia,California,Canary"
 fAtt@CMIP5_model_names = "BCC-CSM1.1,CCSM4,CSIRO-Mk3.6,CanESM2,FGOALS-s2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM,MIROC-ESM-CHEM,MIROC5,MPI-ESM-LR,MRI-CGCM3,NorESM1-M"
 fileattdef(f_output, fAtt )
 SST_regional@unit="K"
 SST_regional!0="years"
 SST_regional!1="regions"
 SST_regional!2="models"
 SST_regional&years = ispan(1900,1999,1)
 SST_regional&regions = ispan(1,5,1)
 SST_regional&models = ispan(1,18,1)
 f_output->SST_regional = SST_regional
