;
; NAME:
;  fsc_NeVIII_06071996.PRO
;
; PURPOSE:
;  To construct full sun 3D data cube from the slit scan of NE VIII

; CALLING SEQUENCE:
;  FULLSUN_NeVIII_06071996

; INPUTS:
;  NONE

; OPTIONAL INPUTS:
;  NONE
;
; OUTPUTS:
;  *sav file with the 3D data cube (also saves to drive)
;  *png. file of Full Disk Peak Intensity of NeVIII [June 7, 1996] (also saves to drive)

; MODIFICATION HISTORY:
;  Written by: Suman Panda 14 Jun 2021 
;  Last Modified: Alyssa Johnson 22 Jun 2021 
; 

;initializes pro file
pro fullsun_NeVIII_06071996

;this keeps track of any files that have dimention other than the 25x360
OpenW,lun,'fulldisk_NeVIII_06071996.txt',/Get_Lun	

;this initializes the array
;the first entry is for NAXIS2
;the second and third entries you get from using CALLING SEQUENCE plot_xy_position_06071996 aka count value
;X and y value of where they are on the Sun
rastor=fltarr(1,3000,3000)
rastor_step=0
ys0=0
ys=0
xs=0
xss=0
strip_no=1

drive='/nfs/hl0/data/ajohnson/CIV_NeVIII_fitsfiles/1996/06/07/'
;fits_files=file_search(drive,'*07704*.fits',count=count)
fits_files=file_search(drive,'*.fits',count=count)
;print,count
;stop	
	
;for i= 0, count-1 do begin
;  xs=xs/3
;  rastor[*,ys:ys+299,xs]=data[0,50:359]

for i=0 ,count-1 do begin ;doing so will give i=0,i=2,i=5
;print,'i=',i
	data= readfits(fits_files(i),Header,/silent)
	
	dimention=size(data)
	if dimention(1) eq 1 and dimention(2) eq 360 then begin ;eq is 1 since NAXIS is 1, NAXIS2 is 360
      		xs0=xss
      		CRVAL2_str=Header(59) ;header column # changes to shave to check with main header at bottom
      		substr1=strsplit(crval2_str,' ',/extract)
      		CRVAL2=FLOAT(substr1(2))
          ;stop
      
      		CRVAL3_str=Header(64) ;header changes so have to check main header at bottom
      		substr2=strsplit(crval3_str,' ',/extract)
      		CRVAL3=FLOAT(substr2(2))
          ;stop
      			
      		xss=ROUND((CRVAL3+1024)) ;1500 gives us the indices for the array (value)
      		yss=ROUND(CRVAL2+1500-60.5) ;ROUND approximates the integer value, ;put (-) the refernce pix  
      		;print,count
          ;stop
          ;try and 
          ;print, xs[300:400]
          ;print, xss[:]
          ;print, ys[:]
          ;print yss[:]
          ;stop
          
          ;print, xs[:]
          ;print, xss[:]
          ;print, ys[:]
          ;print yss[:]
          ;stop	
      		if abs(ys0-yss) ge 250 then begin
      			;print,abs(xs0-xs)
      			;print,xs,xss,xs0,ys,yss
      			ys=yss
      			xs=xss;-abs(xs0-xs)		
      			strip_no=strip_no+1
      			window,xs=1000,ys=1000
      			plot_image,sqrt(max(rastor,dimension=1))
      			cursor,x,y
            ;stop
            ;print,crval2
            ;print,crval3
            ;print,xss
            ;print,yss
            ;find out if data is continuous 
            ;solar x is crval2
            ;solar y is crval3
      		endif				
      		 
          ;print,'xs=',xs
          ;print,'crval2=',crval2
          ;print,'crval3=',crval3  	
      		;print,xs,xss,xs0,ys,yss
      		;print,CRVAL3+1024.,xs,ys	
      		rastor[0,yss:yss+299,xss]=data[0,60:359]
      		ys0=yss
      		
      		if xs ne xss then rastor[0,yss:yss+299,xs]=data[0,60:359]
          ;;;trial ;[0:0] is same as [0]
      
      		if (strip_no MOD 2) eq 1 then xs=xss-1 else xs=xss+1
          ;print,count
          ;stop
	endif else begin
      		printf,lun,fits_files(i)
  endelse



endfor
	window,xs=700,ys=700
	plot_image,sqrt(max(rastor,dimension=1))
  cursor,x,y 


print,xs
print,ys
;stop
fps=30
ovid=idlffvideowrite('sumer_NeVIII_06071996_test.mp4')
xsize=398
ysize=2387
vidstream=ovid.addvideostream(xsize,ysize,fps)
image=tvrd(/true)
;time=ovid.put(vidstream.image)
ovid.cleanup

peak_intensity = max(rastor,dimension=1)
peak_intensity = TRANSPOSE(peak_intensity)
window,xs=1000,ys=1000
;plot_image,sqrt(congrid(peak_intensity,1000,1000))
plot_image,sqrt(max(rastor,dimension=1)),title='Full Disk Peak Intensity in Ne VIII (June 7, 1996)'
write_png,'fullsun_NeVIII_06071996.png',tvrd(/true)





;loadct,0, /silent

;for i =0, xs-1 do begin
;  for j = 0 ny-1 do begin
;      plot_image,transpose(sqrt(max(rastor,dimension=1))),title='Full Disk Peak Intensity in C IV (February 4, 1996)'
;      image=tvrd(/true)
;      time=ovid.put(vidstream,image)
;    endif
;    xs+=1
;  endif
;endfor
;fullsun_NeVIII_06071996 = drive('fsc_NeVIII_06071996.png')
;write_png, 'fsc_NeVIII_06071996.png', TVRD(/TRUE)
;image2 = read_png(fsc_NeVIII_06071996)

;plot_image,sqrt(congrid(peak_intensity[0:2200,300:2500],1000,1000)),title='Full Disk Peak Intensity in Ne VIII (June 7, 1996)',charsize=1.5
;stop



Free_lun,lun

;stop
;fd_Neviii=rastor
;;;Save the 3D data cube as a sav file
;;save,fd_Neviii,filename='fd_Neviii.sav'

;end


;for loop loop initializes and then runs through all fits files
;/silent doesnt print out
;checking for dimension of data and making sure it matches with header (Ex: CIV had larger fits files) ---> good to always check for dimension
;line 83 will print onto file (printf) if any of the data is too large/irrelevant
;/extract is for crvalues (want them in integer values) --->  FLOAT makes integer value
;don't want negative values ---> why sum doing +1024 (line 53)
;all reference pixels ---> we saw that CRVAL2 had 360/2 = 180.5 ---> have to put in correct position
;1 arcsecond in 1 pixel
;using array index (rastor) to make solar y ---> reference pixel changes (check the header to make sure what the reference pixel and what is the difference)



;[SAMPLE HEADER FILE FOR JUNE 7, 1996]

;SIMPLE  =                    T / Written by IDL:  Sat Sep
;BITPIX  =                    8 / Number of bits per data
;NAXIS   =                    2 / Number of data axes ;dimension
;NAXIS1  =                    1 / ;one dimension
;NAXIS2  =                  360 / ;sphere dimension
;EXTEND  =                    T / FITS data may contain ex
;DATE    = '2018-09-08'         /  Creation date of FITS h
;FILENAME= 'sum_19960607_15493794_07704_20.fits' / Name of
;SORIG   = 'GenuineIntel GNU/Linux' / Architecture and OS
;DATASRC = 'Final Distribution (CDROM)' / Data Source
;TELESCOP= 'SOHO    '           /
;INSTRUME= 'SUMER   '           /
;DATE_MOD= '2018-09-08T02:04:01.618' / Last modified
;ORIGIN  = 'SOHO MPS'           / Where Data is Produced
;OBS_SEQ = 'Raster  '           / Name of observing sequen
;DATE_OBS= '1996-06-07T15:49:37.942' / Beginning of Observ
;DATE_END= '1996-06-07T15:49:39.944' / End of Observation
;OBT_TIME=       1212853807.942 / Starting time of acquisi
;OBT_END =       1212853809.944 / End time of acquisition
;LEVEL   = '0       '           / Data processing level
;PRODLVL = 'LZ      '           / Data calibration level
;DETECTOR= 'A       '           / Detector used
;EXPTIME =                 2.00 / Exposure time [s]
;IXWIDTH =              1.00000 / Image width ["]
;IYWIDTH =              300.000 / Image height ["]
;SLIT    = '<2> 1.0" * 300" centered' / SUMER slit
;INSTITUT= '        '           / Name of institutes for c
;CMP_TYPE= 'None    '           / Type of coordinated prog
;CMP_NAME= 'None    '           / Name of campaign observa
;CMP_ID  =                    0 / Campaign number
;STUDY_ID=                  344 / Study number
;STUDY_NM= 'HAS FS 1548/770'    / Study name
;OBJECT  = 'FS      '           / Target
;SCI_OBJ = 'Full sun'           / The science objective
;SCI_SPEC= 'Full Sun in C IV, O IV, Ne VIII HAS' / The spe
;POPUDP  =                    6 / POP/UDP number
;OBS_PROG= '        '           / Name of scientist progra
;SCIENTIS= 'SUMER Team'         / Scientist responsible of
;FFONOFF = 'OFF     '           / On board flat field ON/O
;BINX    =                    1 / Binning X (1 - 1024)
;BINY    =                    3 / Binning Y (1 - 360)
;ROTCMP  =              0.00000 / Solar rotation compensat
;PIMGTYP = 'SINGLE  '           / Processing image Type (S
;COMPRESS=                   16 / Compression method
;COMPAR1 =                    0 / Compression parameter 1
;COMPAR2 =                    0 / Compression parameter 2
;COMPAR3 =                    0 / Compression parameter 3
;DECOMP  =                    T / Decompression (F=No/T=Ye
;FLATCORR=                    F / Flatfield corrected on g
;GEOMCORR=                    F / Geometry corrected on gr
;WCSAXES =                    3 / Number of WCS axis
;CTYPE1  = 'WAVE    '           / Type of CRVAL1
;CUNIT1  = 'Angstrom'           / Axis unit along axis 1
;CRPIX1  =                    1 / Reference pixel along ax
;CRVAL1  =              770.400 / Value at reference pixel
;CDELT1  =            0.0210955 / Axis increments along ax
;CTYPE2  = 'HPLT-TAN'           / Type of CRVAL2
;CUNIT2  = 'arcsec  '           / Axis unit along axis 2
;CRPIX2  =              60.5000 / Reference pixel along ax ;starting point
;CRVAL2  =             -947.100 / Value at reference pixel
;CDELT2  =                    3 / Axis increments along ax
;CTYPE3  = 'HPLN-TAN'           / Type of CRVAL3
;CUNIT3  = 'arcsec  '           / Axis unit along axis 3
;CRPIX3  =                    1 / Reference pixel along ax
;CRVAL3  =             -621.499 / Value at reference pixel
;CDELT3  =              1.00000 / Axis increments along ax
;REFPIX  =                  941 / Original reference Pixel ;org ref pixel IS NOT the ref ref pixel
;DETXSTRT=                  941 / Detector readout start /
;DETXEND =                  941 / Detector readout end / X
;WAVEMIN =              770.400 / Minimum wavelength in im
;WAVEMAX =              770.400 / Maximum wavelength in im
;CORORB  =                    F / Orbitology corrected
;INS_X   =             -624.688 / Pointing of instrument /
;INS_Y   =             -945.000 / Pointing of instrument /
;SOLAR_X =             -621.499 / Instrument pointing / so
;SOLAR_Y =             -947.100 / Instrument pointing / so
;RASTYPE = 'RASTER  '           / Sequence type
;SC_X0   =         -4.78651e-05 / Spacecraft pointing / so
;SC_Y0   =            0.0142032 / Spacecraft pointing / so
;SC_ROLL =             0.193087 / Spacecraft roll angle re
;SOLAR_P0=       -12.7522382029 / Solar angle P0 (degree)
;SOLAR_B0=       0.114591561258 / Solar angle B0 (degree)
;GEIX_OBS=          2.78021e+08 / Geocentric equatorial in
;GEIY_OBS=          1.58241e+09 / Geocentric equatorial in
;GEIZ_OBS=          5.37380e+08 / Geocentric equatorial in
;GSEX_OBS=          1.68582e+09 / Geocentric solar eclipti
;GSEY_OBS=          9.75184e+07 / Geocentric solar eclipti
;GSEZ_OBS=         -1.36429e+08 / Geocentric solar eclipti
;GSMX_OBS=          1.68582e+09 / Geocentric solar magneti
;GSMY_OBS=          1.03438e+08 / Geocentric solar magneti
;GSMZ_OBS=         -1.31997e+08 / Geocentric solar magneti
;HAEX_OBS=         -3.33300e+10 / Heliocentric Aries Eclip
;HAEY_OBS=         -1.46415e+11 / Heliocentric Aries Eclip
;HAEZ_OBS=         -1.43174e+08 / Heliocentric Aries Eclip
;DSUN_OBS=          1.50160e+11 / Distance from Sun
;HGLN_OBS=              0.00000 / Stonyhurst heliographic
;HGLT_OBS=             0.114592 / Stonyhurst heliographic
;CRLN_OBS=              277.426 / Carrington heliographic
;CRLT_OBS=             0.114592 / Carrington heliographic
;XCEN    =             -621.499 / Center of the instrument
;YCEN    =             -947.100 / Center of the instrument
;ANGLE   =             0.193087 / Orientation of instrumen
;UDP_ID  =                  720 / Reference to observing p
;PROG_NM = '        '           / Name of observing progra
;T3TELE  =              53.4400 / Telescope (MC2) temperat
;T3REAR  =              21.6400 / SUMER rear (MC3) tempera
;T3FRONT =              20.1200 / SUMER front (MC4) temper
;T3SPACER=              23.8000 / SUMER spacer (MC6) tempe
;MC2ENC  =                 3197 / SUMER MC2 Encoder Positi
;MC3ENC  =                 2281 / SUMER MC3 Encoder Positi
;MC4ENC  =                 1433 / SUMER MC4 Encoder Positi
;MC6ENC  =                 3383 / SUMER MC6 Encoder Positi
;MC8ENC  =               131072 / SUMER MC8 Encoder Positi
;HKOTIME = '1996-06-07T15:49:34.891' / Time stamp of HK0 R
;
;        /SUMER Raw Image Header Entries
;
;SSYNC0  =                  235 / Header Sync Byte 0 (0xEB
;SSYNC1  =                  144 / Header Sync Byte 1 (0x90
;IMG_REC =                  128 / Header Record Id (0x80)
;SSTYPIMG=                   20 / SUMER image format index
;XSSTAI  =       1212853807.942 / Exposure start time (TAI
;SSOPCNT =                  162 / Operations counter. 0 af
;SSPOPUDP=                    6 / POP/UDP number. 0 no POP
;SSIMGCNT=                    1 / Image counter, 0 at star
;SSLOC   =                   22 / Scientist ID (after Jun
;SSTARGET=                    0 / Target as set by SYS_Ope
;SSFLDATE=                    0 / Never really used
;SSFLREQN=                    0 / Never really used
;SSREFPIX=                  941 / Reference pixel where la
;FREFPIX =                  941 / Ref pixel reduced by det
;SSSTAT  =                   32 / Status
;SSEETRIG=                    0 / Explosive event trigger
;SSFF    =                    0 / Flatfield correction 1=o
;SSFF_T  = 'OFF     '           / Flatfield correction sta
;SSDETTYP=                    1 / Current detector in use
;SSDET_T = 'A       '           / Detector in use as text
;SSINTSTA=                    0 / IIF/TC spectrohelio inte
;SSDETSTA=                  255 / Detector status
;SSSUNY  =                -9995 / SUMER coordinate Y [0.06
;SSSUNZ  =                15120 / SUMER coordinate Z [0.06
;P_SUNY  =             -624.688 / SUMER coordinate Y ["]
;P_SUNZ  =              945.000 / SUMER coordinate Z ["]
;SSEXPTIM=              2.00218 / Exposure time [s]
;SSIIDZ  =                    0 / Inter instrument sun z c
;SSIIDY  =                    0 / Inter instrument sun y c
;SSBPADDY=                   31 / Brightest pix addr Y (sp
;SSBPADDZ=                    6 / Brightest pix addr Z (sp
;SSIMGTOT=                  278 / Total counts in image
;SSROTCOM=              0.00000 / Rotation compensation ti
;SSBPCNTS=                    8 / Brightest pixel counts
;SSACIMGC=                52562 / Accumulative image count
;SSSTEPN =                  414 / Number of raster steps (
;SSSTEPSZ=                   -8 / Raster step size in moto
;SSSLITN =                    2 / Number of slit selected
;SSBINNY =                    1 / Binning factor spectral
;SSBINNZ =                    3 / Binning factor spatial d
;SSXCNT  =                 1884 / X event count (raw value
;SSYCNT  =                 1851 / Y event count (raw value
;X_EVCNT =              13165.4 / X event count (corrected
;Y_EVCNT =              12934.8 / Y event count (corrected
;S_MCPV  =                  238 / High voltage raw value
;S_MCPVF =              5279.00 / High voltage [V]
;SIMCPI  =                  116 / MCP current (Raw value)
;SIMCPIF =              35.7160 / MCP current [uA]
;SSMC2POS=                 4170 / MC2-Azimuth position [0.
;SSMC3POS=                 2540 / MC3-Elev position [0.38"
;SSMC4POS=                 -954 / MC4-Slit select position
;SSMC6POS=                18664 / MC6-Grat position [half
;SSMC8POS=                10054 / MC8-Scan mirror position
;SSMCERR =                    0 / Motor controller error
;SSCOMPRM=                   16 / Compression method
;SSWAVEL =              770.400 / Wavelength at reference
;SSCOMPP1=                    0 / Compression parameter 1
;SSCOMPP2=                    0 / Compression parameter 2
;SSCOMPP3=                    0 / Compression parameter 3
;SSADMCNT=                  720 / Admin cnt (gives the dat
;XSSMDU  =                    0 / Missing data in image 0=
;XSSQAC  =                    0 / Quality of image data 0=
;XSSFID  =                    0 / File id from TM file cat
;XSSFPTR =             54439680 / Pointer to image positio
;XSCDID  =                   60 / CD ID of TM file
;XSSEQID =                    5 / CD Sequence ID of TM fil
;XSTMFILE=             61591204 / TM filename without ext
;PROCVERS=                  552 / SVN Version of write_sum
;
;        / Data integrity check block
;
;CHECKSUM= '3Ca4A9U25AZ2A9Z2'   / HDU checksum updated 201
;DATASUM = '875967540'          / data unit checksum updat
;COMMENT FITS (Flexible Image Transport System) format is
;COMMENT and Astrophysics', volume 376, page 359; bibcode
;COMMENT The original raw image header is stored in the ex
;COMMENT The wavelength is stored low -> high
;COMMENT The data is stored with north up!
;COMMENT For documentation on SUMER FITS Data see:
;COMMENT Solar Physics June 2014, Volume 289, Issue 6, pp
;HISTORY Written with write_sumerfits Version:$Revision: 5
;HISTORY Applied decomp_method16 $Revision: 244 $
;HISTORY Applied update_fits $Revision: 561 $  8-Sep-2018
END
