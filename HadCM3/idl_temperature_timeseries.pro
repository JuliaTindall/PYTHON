pro temperature_timeseries_frompg

;have started amending this to add in temperatures from HadCM3
;hope I have not mucked it up.



; this program will plot the timeseries from the pg files to
; investigate the drift
;
close,/all


expt='xhcph'
HadCM3='y'
database='y'


set_plot,'ps'
filename=STRING('/nfs/hera1/earjcti/IDLPLOTS/HadGEM/',expt,'_temperature_drift.ps')
filename=STRCOMPRESS(filename,/remove_all)
device,filename=filename,/landscape
!p.multi=[0,1,2]

;yearstart=501
;yearend=600
yearstart=01
yearend=499

centstart='u'

; find place in array for centstart
allcents=['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0']
ncents=n_elements(allcents)
for cent=0,ncents-1 do begin
   if allcents(cent) eq centstart then begin
      centstart_ss=cent
   endif
endfor



tsize=yearend-yearstart+1
print,'tsize is',tsize
if HadCM3 eq 'y' then begin
   xsize=288
   ysize=144
   zsize=20
endif else begin
   xsize=360
   ysize=216
   zsize=40
endelse

if database eq 'n' then begin
   alldata=fltarr(xsize,ysize,zsize,tsize)
   globaltemp=fltarr(zsize,tsize)
   totaltemp=fltarr(tsize) 
   correct=intarr(tsize)

  
   for t=yearstart,yearend do begin 
      century=FLOOR(t/100.)
      extrause=allcents(centstart_ss+century)
      yearuse=t-(century * 100)
      
      if yearuse lt 10 then begin
         yearuse=STRING('0',yearuse)
         yearuse=STRCOMPRESS(yearuse,/remove_all)
      endif
      
      
      plotdata,expt,alldata,yearuse,extrause,yearuse,depths,$
               glob_trac,yearstart,latitudes,depthdiff,fulltemp,requires_correction,$
               HadCM3,xsize,ysize,zsize
      
      globaltemp(*,t-yearstart)=glob_trac(*)
      totaltemp(t-yearstart)=fulltemp
      if requires_correction eq 'y' then begin
         correct(t-yearstart)=1
      endif
      
   endfor
endif else begin ; use experiment from database
   filename='/nfs/hera2/apps/metadata/experiments/'+expt+'/timeseries/'+expt+'o@pg_OceanTemperature.nc'

   ncdfid=ncdf_open(filename,/nowrite)
; get dimensions
   latid=ncdf_varid(ncdfid,'latitude')
   longid=ncdf_varid(ncdfid,'longitude')
   depthid=ncdf_varid(ncdfid,'depth')
   timeid=ncdf_varid(ncdfid,'t')

   ncdf_varget,ncdfid,latid,latitudes
   ncdf_varget,ncdfid,longid,longitudes
   ncdf_varget,ncdfid,depthid,depth
   ncdf_varget,ncdfid,timeid,times

   xsize=N_ELEMENTS(longitudes)
   ysize=N_ELEMENTS(latitudes)
   zsize=N_ELEMENTS(depth)
   tsize=N_ELEMENTS(times)

   alldata=fltarr(xsize,ysize,zsize,tsize)
   globaltemp=fltarr(zsize,tsize)
   totaltemp=fltarr(tsize)


    varid=ncdf_varid(ncdfid,'insitu_T')
    ncdf_varget,ncdfid,varid,alldata
    ncdf_close,ncdfid
    
    ; calculate depth difference
    depthdiff=fltarr(zsize)
    depthdiff(0)=depth(0)*2.
    for i=1,zsize-2 do begin
       depthdiff(i)=0.5 * (depth(i+1)-depth(i-1))
    endfor
    depthdiff(zsize-1)=(depth(zsize-1)-depth(zsize-2))*2.

    ; get global ratio for each time
    globaltemp=fltarr(zsize,tsize)
    totT=fltarr(zsize,tsize)
    div=fltarr(zsize,tsize)
    fulltemp=fltarr(tsize)
    fulldiv=fltarr(tsize)
    for t=0,tsize-1 do begin
       print,'t is',t
       for k=0,zsize-1 do begin
          for j=0,ysize-1 do begin
             for i=0,xsize-1 do begin
                if alldata(i,j,k,t) lt 1.0E10 then begin
                   totT(k,t)=totT(k,t)+(alldata(i,j,k,t)*cos(latitudes(j) * 2. * !pi/360.))
                   div(k,t)=div(k,t)+cos(latitudes(j) * 2. * !pi/360.)
                   fulltemp(t)=fulltemp(t)+(alldata(i,j,k,t)*cos(latitudes(j) * 2. * !pi/360.)*depthdiff(k))
                   fulldiv(t)=fulldiv(t)+(cos(latitudes(j) * 2. * !pi/360.)*depthdiff(k))
                endif
             endfor
          endfor
       endfor
    endfor

      
    globaltemp=totT/div
    totaltemp=fulltemp/fulldiv
    correct=intarr(tsize) ; don't correct unless its clear we need to 

    
endelse
 
; correct if required
correct_diff=fltarr(zsize)
if correct(0) eq 1 then begin
print,'we are correcting'
    for t=1,tsize-1 do begin
        if correct(t) ne 1 then begin 
            correct_at=t
            for k=0,zsize-1 do begin
                correct_diff(k)=globaltemp(k,t)-globaltemp(k,t-1)
            endfor
            correct_all=totaltemp(t)-totaltemp(t-1)
print,correct_diff(0),globaltemp(0,t),globaltemp(0,t-1)
print,correct_all,totaltemp(t),totaltemp(t-1)
            goto,foundval
        endif
    endfor
    foundval:
    for t=0,tsize-1 do begin
        if correct(t) eq 1 then begin
            globaltemp(*,t)=globaltemp(*,t)+correct_diff(*)
            totaltemp(t)=totaltemp(t)+correct_all
        endif
        if correct(t) eq 0 then begin
            goto,corrected
        endif
    endfor
    corrected:
endif
            



for k=0,zsize-1 do begin
    titlename=STRING('global temperature for level',k)
    plot,indgen(tsize),globaltemp(k,*),title=titlename,xtitle=year,ytitle='degC',ystyle=1
endfor



plot,indgen(tsize),totaltemp,title='TOTAL OCEAN TEMPERATURE',xtitle=year,ytitle='degC',ystyle=1

textout='/nfs/hera1/earjcti/IDLPLOTS/HadGEM/textfiles/'+expt+'drift.tex'
openw,1,textout

printf,1,'file produced by temperature_timeseries_frompg'
printf,1,'total drift'
for i=0,tsize-1 do begin
    printf,1,i+yearstart,totaltemp(i)
endfor
printf,1,'levels 0-4'
for i=0,tsize-1 do begin
    printf,1,i+yearstart,globaltemp(0,i),globaltemp(1,i),globaltemp(2,i),globaltemp(3,i),globaltemp(4,i)
endfor
printf,1,'levels 5-9'
for i=0,tsize-1 do begin
    printf,1,i+yearstart,globaltemp(5,i),globaltemp(6,i),globaltemp(7,i),globaltemp(8,i),globaltemp(9,i)
endfor
printf,1,'levels 10-14'
for i=0,tsize-1 do begin
    printf,1,i,globaltemp(10,i),globaltemp(11,i),globaltemp(12,i),globaltemp(13,i),globaltemp(14,i)
endfor
printf,1,'levels 15-19'
for i=0,tsize-1 do begin
    printf,1,i+yearstart,globaltemp(15,i),globaltemp(16,i),globaltemp(17,i),globaltemp(18,i),globaltemp(19,i)
endfor
if HadCM3 ne 'y' then begin
   printf,1,'levels 20-24'
   for i=0,tsize-1 do begin
      printf,1,i+yearstart,globaltemp(20,i),globaltemp(21,i),globaltemp(22,i),globaltemp(23,i),globaltemp(24,i)
   endfor
   printf,1,'levels 25-29'
   for i=0,tsize-1 do begin
      printf,1,i+yearstart,globaltemp(25,i),globaltemp(26,i),globaltemp(27,i),globaltemp(28,i),globaltemp(29,i)
   endfor
   printf,1,'levels 30-34'
   for i=0,tsize-1 do begin
      printf,1,i+yearstart,globaltemp(30,i),globaltemp(31,i),globaltemp(32,i),globaltemp(33,i),globaltemp(34,i)
   endfor
   printf,1,'levels 35-39'
   for i=0,tsize-1 do begin
      printf,1,i+yearstart,globaltemp(35,i),globaltemp(36,i),globaltemp(37,i),globaltemp(38,i),globaltemp(39,i)
   endfor
endif

close,1



;oplot,indgen(tsize),smoothval


device,/close

print,'end of program'
end

;============================================
pro plotdata,expt,alldata,yearuse,extrause,t,depth,globaltemp,$
             yearstart,latitudes,depthdiff,fulltemp,requires_correction,$
             HadCM3,xsize,ysize,zsize

;

; 1. get data


if HadCM3 eq 'y' then begin
   filein=STRING('/nfs/hera1/earjcti/um/netcdf/',expt,'_netcdf/',expt,'o@pg',extrause,yearuse,'c1.nc')
endif else begin
   filein=STRING('/nfs/hera1/earjcti/um/HadGEM_data/',expt,'/netcdf/',expt,'o@pg',extrause,yearuse,'c1.nc')
endelse
filein=STRCOMPRESS(filein,/remove_all)
print,filein

; check if file exists if it doesn't use 'l' version'
a=findfile(filein,count=count)
print,'count is',count	
requires_correction='n'
if count eq 0 then begin
    filein=STRING('/nfs/hera1/earjcti/um/',expt,'/netcdf/',expt,'o@da',extrause,yearuse,'c1.nc')
    filein=STRCOMPRESS(filein,/remove_all)
    requires_correction='y'  ; if we are using @da files we need to add a 
                                ; correction factor as we are
                                ; comparing a single point from
                                ; december with the whole year
endif 

print,filein,extrause


ncdfid=ncdf_open(filein,/nowrite)
; get dimensions
latid=ncdf_varid(ncdfid,'latitude')
if count eq 0 then begin
    longid=ncdf_varid(ncdfid,'longitude')
    depthid=ncdf_varid(ncdfid,'depth')
endif else begin
    longid=ncdf_varid(ncdfid,'longitude')
    depthid=ncdf_varid(ncdfid,'depth_1')
endelse
ncdf_varget,ncdfid,latid,latitudes
ncdf_varget,ncdfid,longid,longitudes
ncdf_varget,ncdfid,depthid,depth
xsize=N_ELEMENTS(longitudes)
ysize=N_ELEMENTS(latitudes)
zsize=N_ELEMENTS(depth)

if count eq 0 then begin    
    varid=ncdf_varid(ncdfid,'temp')
endif else begin
    varid=ncdf_varid(ncdfid,'temp')
endelse
ncdf_varget,ncdfid,varid,temp
ncdf_close,ncdfid
print,t,t-yearstart,N_ELEMENTS(alldata(0,0,*))

if count eq 0 and xsize eq 362 then xsize=360

for i=0,xsize-1 do begin
    for j=0,ysize-1 do begin
       for k=0,zsize-1 do begin
           alldata(i,j,k,t-yearstart)=temp(i,j,k,0)
       endfor
    endfor
endfor
    
; calculate depth difference
depthdiff=fltarr(zsize)
depthdiff(0)=depth(0)*2.
for i=1,zsize-2 do begin
	depthdiff(i)=0.5 * (depth(i+1)-depth(i-1))
endfor
depthdiff(zsize-1)=(depth(zsize-1)-depth(zsize-2))*2.



; get global ratio for each time
globaltemp=fltarr(zsize)
totT=fltarr(zsize)
div=fltarr(zsize)
fulltemp=0
fulldiv=0
for j=0,ysize-1 do begin
    for i=0,xsize-1 do begin
        for k=0,zsize-1 do begin
            if temp(i,j,k,0) lt 1.0E10 then begin
                totT(k)=totT(k)+(temp(i,j,k,0)*cos(latitudes(j) * 2. * !pi/360.))
                div(k)=div(k)+cos(latitudes(j) * 2. * !pi/360.)
                fulltemp=fulltemp+(temp(i,j,k,0)*cos(latitudes(j) * 2. * !pi/360.)*depthdiff(k))
                fulldiv=fulldiv+(cos(latitudes(j) * 2. * !pi/360.)*depthdiff(k))
            endif
        endfor
    endfor
endfor
globaltemp=totT/div
fulltemp=fulltemp/fulldiv



return


end
