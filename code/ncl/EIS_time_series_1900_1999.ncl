; This ncl program computes annual-mean EIS time series based on 20th-century simulations with 20 CMIP5 models
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

 Lv  = 2470000. ; J kg-1
 Cp  =  1015.   ; J kg K-1
 gr  =  9.8     ; m/s2
 Rv  = 461.5    ; J/kg K
 Ra  = 286.9    ; J/kg K
 EIS_regional = new((/100, 5,18/),double) 

 diri_init  = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo")
 model_index = (/2,5,7,8,10,11,12,13,15,16,17,19,22,23,24,25,28,29/)
 diri = diri_init(model_index)

 start_year   = (/1850,1850,1850,1850,1850,1860,1861,1896,1850,1885,1885,1850,1850,1850,1850,1850,1850,1850/)
 length1 = 1999-1900+1

 do model_no = 0,17;
; historical
  file_ta_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/hist/ta_*r1i1p1*.nc")
  file_ps_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/hist/ps_*.nc")
  file_tas_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/hist/tas_*.nc")

  f_ta_hist_0    = addfile(file_ta_hist(0), "r")   ; note the "s" of addfile
  f_ta_hist      = addfiles(file_ta_hist, "r")   ; note the "s" of addfile
  f_tas_hist     = addfiles(file_tas_hist, "r")   ; note the "s" of addfile
  f_ps_hist      = addfiles(file_ps_hist, "r")   ; note the "s" of addfile
  lat  = f_ta_hist_0->lat
  lon  = f_ta_hist_0->lon
  lev  = f_ta_hist_0->plev
  EIS_3d = new((/12,100,dimsizes(lat),dimsizes(lon)/),double)

  do i=0,length1-1
    print ("model_no="+model_no +"    year="+i)
    do imonth=0,11
        aa=minind(abs(lev-70000));
        ta_djfmam_mean = f_ta_hist[:]->ta((1900-start_year(model_no)+i)*12+imonth,aa,:,:);
        tas_djfmam_mean = f_tas_hist[:]->tas((1900-start_year(model_no)+i)*12+imonth,:,:);
        sp = f_ps_hist[:]->ps((1900-start_year(model_no)+i)*12+imonth,:,:);
        if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES") ; spatial interpolation is needed
            f_tas_hist_0    = addfile(file_tas_hist(0), "r")   ; note the "s" of addfile
            tas_lat         = f_tas_hist_0->lat
            tas_lon         = f_tas_hist_0->lon
            tas_new         = linint2(tas_lon,tas_lat,tas_djfmam_mean,True,lon,lat,0)
            delete(tas_djfmam_mean)
            tas_djfmam_mean = tas_new
            delete(tas_new)
            sp_new          = linint2(tas_lon,tas_lat,sp,True,lon,lat,0)
            delete(sp)
            sp = sp_new
            delete(sp_new)
        end if
        LTS = ta_djfmam_mean(:,:)*((1000./700.)^(Ra/Cp)) - tas_djfmam_mean(:,:)*((100000./sp(:,:))^(Ra/Cp)) 
        T850 = (ta_djfmam_mean(:,:)+tas_djfmam_mean(:,:))/2.0 - 273.15 
        tas  = tas_djfmam_mean(:,:)
        es   = 6.11*exp(17.3*T850/(T850+237.3))     ; surface saturation vapor pressure
        qs   = 0.622*es/(850-es)            ; mixing ratio in surface
        Gamma = (gr/Cp)*(1-(1+Lv*qs/Ra/(T850+273.15))/(1+Lv*Lv*qs/Cp/Rv/(T850+273.15)/(T850+273.15)))
        z700  = Ra*tas_djfmam_mean(:,:)*log(sp(:,:)/70000.)/gr
        rh    = tas_djfmam_mean(:,:)
        rh(:,:)  = 80 
        tdk   = dewtemp_trh(tas_djfmam_mean(:,:),rh)
        plcl  = lclvl(sp(:,:),tas_djfmam_mean(:,:),tdk)
        zlcl  = Ra*tas_djfmam_mean(:,:)*log(sp(:,:)/plcl)/gr
        EIS   = LTS - Gamma*(z700-zlcl)
        EIS_3d(imonth,i,:,:) = EIS 
        delete(ta_djfmam_mean)
        delete(tas_djfmam_mean)
        delete(sp)
        delete(es)
        delete(qs)
        delete(Gamma)
        delete(z700)
        delete(rh)
        delete(tdk)
        delete(plcl)
        delete(zlcl)
        delete(LTS)
        delete(T850)
        delete(tas)
        delete(EIS)
   end do
  end do;
  delete(file_tas_hist)
  delete(file_ta_hist)
  delete(file_ps_hist)
  delete(f_tas_hist)
  delete(f_ps_hist)
  delete(f_ta_hist)
  delete(f_ta_hist_0)
  if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
     delete(f_tas_hist_0)
  end if

  file_sftlf = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/sftlf*")
  f_cl_sftlf_0   = addfile(file_sftlf(0), "r")   ; note the "s" of addfile
  sftlf = f_cl_sftlf_0->sftlf
  if (diri(model_no).eq."MIROC-ESM" .or. diri(model_no).eq."MIROC-ESM-CHEM")
     sftlf  = sftlf*100
  end if
  if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
     sftlf_new         = linint2(tas_lon,tas_lat,sftlf,True,lon,lat,0)
     delete(sftlf)
     sftlf = sftlf_new
     delete(sftlf_new)
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
     EIS0   = dim_avg_n(EIS_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lon_2d.ge.360-110 .and. lon_2d.le.360-70 .and. lat_2d.ge.-30 .and. lat_2d.le.-10 .and. sftlf.lt.0.01, EIS0,-999)
     weight = where(lon_2d.ge.360-110 .and. lon_2d.le.360-70 .and. lat_2d.ge.-30 .and. lat_2d.le.-10 .and. sftlf.lt.0.01, weight, -999)
     weight@_FillValue = -999
     EIS0@_FillValue = -999
     EIS_regional(iyear,0,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)

; Namibia
     EIS0   = dim_avg_n(EIS_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .ge. 360-25. .or. lon_2d.le.15) .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .ge. 360-25. .or. lon_2d.le.15) .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     EIS_regional(iyear,1,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)

; Australia
     EIS0   = dim_avg_n(EIS_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     EIS_regional(iyear,2,model_no) = sum(EIS0*weight)/sum(weight)  
     delete(EIS0)
     delete(weight); 

; California
     EIS0   = dim_avg_n(EIS_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and.  sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     EIS_regional(iyear,3,model_no) = sum(EIS0*weight)/sum(weight)  
     delete(EIS0)
     delete(weight)

; Canary
     EIS0   = dim_avg_n(EIS_3d(:,iyear,:,:),0)
     weight = cos(3.14*lat_2d/180.)
     EIS0 = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, EIS0, -999)
     weight = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, weight, -999)
     EIS0@_FillValue = -999
     weight@_FillValue = -999
     EIS_regional(iyear,4,model_no) = sum(EIS0*weight)/sum(weight)
     delete(EIS0)
     delete(weight)
  end do
  delete(sftlf)
  delete(lat)
  delete(lon)
  if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
    delete(tas_lat)
    delete(tas_lon)
  end if
  delete(lev)
  delete(lat_2d)
  delete(lon_2d)
  delete(aa)
  delete(EIS_3d)
 end do;

 f_output = addfile("EIS_time_series_1900_1999_CMIP5.nc","w");  
 fAtt               = True            ; assign file attributes
 fAtt@title         = "annual-mean EIS time series in 18 CMIP5 models during the period of 1900-1999"
 fAtt@creation_date = systemfunc ("date")
 fAtt@region_names = "Peru,Namibia,Australia,California,Canary"
 fAtt@CMIP5_model_names = "BCC-CSM1.1,CCSM4,CSIRO-Mk3.6,CanESM2,FGOALS-s2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM,MIROC-ESM-CHEM,MIROC5,MPI-ESM-LR,MRI-CGCM3,NorESM1-M"
 fileattdef(f_output, fAtt )
 EIS_regional@unit="K"
 EIS_regional!0="years"
 EIS_regional!1="regions"
 EIS_regional!2="models"
 EIS_regional&years = ispan(1900,1999,1)
 EIS_regional&regions = ispan(1,5,1)
 EIS_regional&models = ispan(1,18,1)
 f_output->EIS_regional = EIS_regional


