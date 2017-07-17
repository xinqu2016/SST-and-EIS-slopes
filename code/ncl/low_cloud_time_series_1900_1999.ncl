; This ncl program computes annual-mean time series of low cloud cover based on 20th-century simulations with 20 CMIP5 models
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; data directory
 diri_init  = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo")
 model_index = (/2,5,7,8,10,11,12,13,15,16,17,19,22,23,24,25,28,29/)
 diri = diri_init(model_index)
 start_year   = (/1850,1850,1850,1850,1850,1860,1861,1896,1850,1885,1885,1850,1850,1850,1850,1850,1850,1850/)
 length1 = 1999-1900+1

 low_cl_regional = new((/100,5,18/),double) ; 100 years, 5 regions and 18 models

 do model_no = 0,17
; historical
  file_cl_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/hist/cl_*")
  f_cl_hist_0    = addfile(file_cl_hist(0), "r")   ; note the "s" of addfile
  f_cl_hist      = addfiles(file_cl_hist, "r")   ; note the "s" of addfile
  lat  = f_cl_hist_0->lat
  lon  = f_cl_hist_0->lon
  b    = f_cl_hist_0->b

  lcc = new((/12,100,dimsizes(lat),dimsizes(lon)/),double) 
  lcc(:,:,:,:)=0
  if (diri(model_no).eq. "CanESM2" .or. diri(model_no).eq."IPSL-CM5A-LR" .or. diri(model_no).eq."MPI-ESM-LR")
    a    = f_cl_hist_0->ap
  else 
;     HadGEM2 uses z-coordinate so we convert it into p-coordinate
     if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
        a    = f_cl_hist_0->lev
        orog = f_cl_hist_0->orog
        file_tas_hist = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/"+diri(model_no)+"/hist/tas_*nc")
        f_tas_hist    = addfiles(file_tas_hist,"r")
        file_ps_hist  = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/"+diri(model_no)+"/hist/ps_*nc")
        f_ps_hist     = addfiles(file_ps_hist,"r")
     else
        a    = f_cl_hist_0->a*f_cl_hist_0->p0
        if(diri(model_no).eq."FGOALS-s2")
             a  = f_cl_hist_0->a*f_cl_hist_0->p0/1e05
        end if
     end if
   end if

  do i=0,length1-1
    print ("model_no="+model_no +"    year="+i)
    do imonth=0,11
        cl_djfmam_mean = f_cl_hist[:]->cl((1900-start_year(model_no)+i)*12+imonth,:,:,:);
        if (diri(model_no).eq."CCSM4")
           cl_djfmam_mean = cl_djfmam_mean*100
        end if
        if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
            tas_djfmam_mean = f_tas_hist[:]->tas((1900-start_year(model_no)+i)*12+imonth,:,:);
            ps_djfmam_mean  = f_ps_hist[:]->ps((1900-start_year(model_no)+i)*12+imonth,:,:);
        else
            ps_djfmam_mean = f_cl_hist[:]->ps((1900-start_year(model_no)+i)*12+imonth,:,:);
        end if

        pr_3d = new((/dimsizes(b),dimsizes(lat),dimsizes(lon)/),double)
        if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
               do ilev=0, dimsizes(b)-1
                     pr_3d(ilev,:,:) = a(ilev)+b(ilev)*orog(:,:)
                     pr_3d(ilev,:,:) = ps_djfmam_mean(:,:)*exp(-(pr_3d(ilev,:,:)*9.8)/(287.*tas_djfmam_mean(:,:)))
               end do
               delete(tas_djfmam_mean)
        else
               do ilev=0, dimsizes(b)-1
                     pr_3d(ilev,:,:) =a(ilev)+b(ilev)*ps_djfmam_mean(:,:)
               end do
        end if
        delete(ps_djfmam_mean)
        cl_djfmam_mean = where(pr_3d.ge.68000,cl_djfmam_mean/100.,0)
        cc_temp = cl_djfmam_mean(0,:,:);
        do ilev=1,dimsizes(b)-1 ; random overlap assumption 
            cc_temp(:,:) = cc_temp(:,:)+cl_djfmam_mean(ilev,:,:)-cc_temp(:,:)*cl_djfmam_mean(ilev,:,:);
        end do;
        lcc(imonth,i,:,:)= 100*cc_temp;
        delete(cc_temp)
        delete(cl_djfmam_mean)
        delete(pr_3d)
   end do
  end do
  delete(file_cl_hist)
  delete(f_cl_hist)
  delete(f_cl_hist_0)
  if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
    delete(file_tas_hist)
    delete(file_ps_hist)
    delete(f_tas_hist)
    delete(f_ps_hist)
  end if 

 file_sftlf = systemfunc("ls /Volumes/Data_Server/CMIP5_modeloutput/mo/" + diri(model_no) + "/sftlf*")
 f_cl_sftlf_0   = addfile(file_sftlf(0), "r")   ; note the "s" of addfile
 sftlf = f_cl_sftlf_0->sftlf
 if (diri(model_no).eq."MIROC-ESM" .or. diri(model_no).eq."MIROC-ESM-CHEM")
     sftlf  = sftlf*100
 end if

 lat_2d  = new((/dimsizes(lat),dimsizes(lon)/),double)
 lon_2d  = new((/dimsizes(lat),dimsizes(lon)/),double)
 weight_2d  = new((/dimsizes(lat),dimsizes(lon)/),double)
 do k=0,dimsizes(lon)-1
     lat_2d(:,k)=lat(:)
 end do
 weight_2d = cos(3.1415*lat_2d/180.)
 do k=0,dimsizes(lat)-1
     lon_2d(k,:)=lon(:)
 end do
 do i=0,length1-1
     low_cl_djfmam_mean = dim_avg_n(lcc(:,i,:,:),0)
;    Peru
     low_cl_djfmam_mean_0 = where(lat_2d .le.-10. .and. lat_2d .ge.-30. .and. lon_2d .le. 360-70. .and. lon_2d .ge. 360.-110 .and. sftlf.le.0.01, low_cl_djfmam_mean, -999)
     weight_2d_0          = where(lat_2d .le.-10. .and. lat_2d .ge.-30. .and. lon_2d .le. 360-70. .and. lon_2d .ge. 360.-110 .and. sftlf.le.0.01, weight_2d, -999)
     low_cl_djfmam_mean_0@_FillValue = -999
     weight_2d_0@_FillValue = -999
     low_cl_regional(i,0,model_no) = sum(low_cl_djfmam_mean_0*weight_2d_0)/sum(weight_2d_0)
     delete(low_cl_djfmam_mean_0)
     delete(weight_2d_0)

;    Namibia
     low_cl_djfmam_mean_0 = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .le. 15. .or. lon_2d.ge.360-25) .and. sftlf.le.0.01, low_cl_djfmam_mean, -999)
     weight_2d_0 = where(lat_2d .le.-7.5 .and. lat_2d .ge.-27.5 .and. (lon_2d .le. 15. .or. lon_2d.ge.360-25) .and. sftlf.le.0.01, weight_2d, -999)
     low_cl_djfmam_mean_0@_FillValue = -999
     weight_2d_0@_FillValue = -999
     low_cl_regional(i,1,model_no) = sum(low_cl_djfmam_mean_0*weight_2d_0)/sum(weight_2d_0)
     delete(low_cl_djfmam_mean_0)
     delete(weight_2d_0)

;    Australia
     low_cl_djfmam_mean_0 = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, low_cl_djfmam_mean, -999)
     weight_2d_0 = where(lat_2d .le.-17.5 .and. lat_2d .ge.-37.5 .and. lon_2d .ge. 75. .and. lon_2d .le. 115 .and. sftlf.le.0.01, weight_2d, -999)
     low_cl_djfmam_mean_0@_FillValue = -999
     weight_2d_0@_FillValue = -999
     low_cl_regional(i,2,model_no) = sum(low_cl_djfmam_mean_0*weight_2d_0)/sum(weight_2d_0)
     delete(low_cl_djfmam_mean_0)
     delete(weight_2d_0)

;    California
     low_cl_djfmam_mean_0 = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and. sftlf.le.0.01, low_cl_djfmam_mean, -999)
     weight_2d_0 = where(lat_2d .le.35 .and. lat_2d .ge.15 .and. lon_2d .ge. 360-155 .and. lon_2d .le. 360-115. .and. sftlf.le.0.01, weight_2d, -999)
     low_cl_djfmam_mean_0@_FillValue = -999
     weight_2d_0@_FillValue = -999
     low_cl_regional(i,3,model_no) = sum(low_cl_djfmam_mean_0*weight_2d_0)/sum(weight_2d_0)
     delete(low_cl_djfmam_mean_0)
     delete(weight_2d_0)

;    Canary
     low_cl_djfmam_mean_0 = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, low_cl_djfmam_mean, -999)
     weight_2d_0 = where(lat_2d .le.30 .and. lat_2d .ge.10 .and. lon_2d .ge. 360-55 .and. lon_2d .le. 360-15 .and. sftlf.le.0.01, weight_2d, -999)
     low_cl_djfmam_mean_0@_FillValue = -999
     weight_2d_0@_FillValue = -999
     low_cl_regional(i,4,model_no) = sum(low_cl_djfmam_mean_0*weight_2d_0)/sum(weight_2d_0)
     delete(low_cl_djfmam_mean_0)
     delete(weight_2d_0)
     delete(low_cl_djfmam_mean)
 end do

 delete(lat)
 delete(lon)
 delete(lat_2d)
 delete(lon_2d)
 delete(weight_2d)
 delete(lcc)
 delete(a)
 delete(b)
 delete(sftlf)
 if (diri(model_no).eq."HadGEM2-CC" .or. diri(model_no).eq."HadGEM2-ES")
    delete(orog)
  end if 
end do 
f_output = addfile("low_cloud_time_series_1900_1999_CMIP5.nc","c");
fAtt               = True            ; assign file attributes
fAtt@title         = "annual-mean LCC time series in 18 CMIP5 models during the period of 1900-1999"
fAtt@creation_date = systemfunc ("date")
fAtt@region_names = "Peru,Namibia,Australia,California,Canary"
fAtt@CMIP5_model_names = "BCC-CSM1.1,CCSM4,CSIRO-Mk3.6,CanESM2,FGOALS-s2,GFDL-CM3,GFDL-ESM2G,GFDL-ESM2M,GISS-E2-R,HadGEM2-CC,HadGEM2-ES,IPSL-CM5A-LR,MIROC-ESM,MIROC-ESM-CHEM,MIROC5,MPI-ESM-LR,MRI-CGCM3,NorESM1-M"
fileattdef(f_output, fAtt )
low_cl_regional@unit="%"
low_cl_regional!0="years"
low_cl_regional!1="regions"
low_cl_regional!2="models"
low_cl_regional&years = ispan(1900,1999,1)
low_cl_regional&regions = ispan(1,5,1)
low_cl_regional&models = ispan(1,18,1)
f_output->low_cl_regional = low_cl_regional

