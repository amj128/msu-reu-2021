; NAME:
;  fsc_NeVIII_02021996.PRO
;
; PURPOSE:
;  To construct full sun 3D data cube from the slit scan of NE VIII

; CALLING SEQUENCE:
;  FULLSUN_NeVIII_02021996

; INPUTS:
;  NONE

; OPTIONAL INPUTS:
;  NONE
;
; OUTPUTS:
;  *sav file with the 3D data cube (also saves to drive)
;  *png. file of Full Disk Peak Intensity of NeVIII [February 2, 1996] (also saves to drive)

; MODIFICATION HISTORY:
;  Written by: Suman Panda 14 Jun 2021 
;  Last Modified: Alyssa Johnson 16 Jun 2021 
; 


;initializes pro file
pro fullsun_02021996

;this keeps track of any files that have dimention other than the 2x360
OpenW,lun,'fulldisk_NeVIII_02021996.txt',/Get_Lun	

;this initializes the array
;the first entry is for NAXIS2
;the second and third entries you get from using CALLING SEQUENCE plot_xy_position_02021996 aka count value
rastor=fltarr(2,7500,7500)
rastor_step=0
ys0=0
ys=0
xs=0
xss=0
strip_no=1

drive='/disk/data/ajohnson/CIV_NeVIII_fitsfiles/1996/02/02/'
fits_files=file_search(drive,'*.fits',count=count)	
;stop	
	
for i= 0, count-1 do begin
	data= readfits(fits_files(i),Header,/silent)
	
	dimention=size(data)
	if dimention(1) eq 25 and dimention(2) eq 360 then begin
		xs0=xss
		CRVAL2_str=Header(68)
		substr1=strsplit(crval2_str,' ',/extract)
		CRVAL2=FLOAT(substr1(2))

		CRVAL3_str=Header(73)
		substr2=strsplit(crval3_str,' ',/extract)
		CRVAL3=FLOAT(substr2(2))
			
		xss=ROUND((CRVAL3+1024))
		yss=ROUND(CRVAL2+1500-830)
		
		if abs(ys0-yss) ge 250 then begin
			print,abs(xs0-xs)
			;print,xs,xss,xs0,ys,yss
			ys=yss
			xs=xss;-abs(xs0-xs)		
			strip_no=strip_no+1
			window,xs=1000,ys=1000
			plot_image,sqrt(congrid(max(rastor,dimension=1),1000,1000))
			cursor,x,y	
		endif				
		 	
		print,xs,xss,xs0,ys,yss
		print,CRVAL3+1024.,xs,ys	
		rastor[0:24,yss:yss+299,xss]= data[0:24,60:359]
		ys0=yss
		
		if xs ne xss then rastor[0:24,yss:yss+299,xs]=data[0:24,60:359];;;trial

		if (strip_no MOD 2) eq 1 then xs=xss-1 else xs=xss+1
	endif else begin
		printf,lun,fits_files(i)
	endelse
	
endfor
	window,xs=1000,ys=1000
	tvscl,sqrt(congrid(max(rastor,dimension=1),1000,1000))
	cursor,x,y
;stop

peak_intensity = max(rastor,dimension=1)
peak_intensity = TRANSPOSE(peak_intensity)
window,xs=1000,ys=1000
plot_image,sqrt(congrid(peak_intensity,1000,1000))

;fullsun_NeVIII_06071996 = drive('fullsun_NeVIII_02021996.png')
;write_png, 'fullsun_NeVIII_02021996.png', TVRD(/TRUE)
;image2 = read_png(fullsun_NeVIII_02021996)

plot_image,sqrt(congrid(peak_intensity[0:2200,300:2500],1000,1000)),title='Full Disk Peak Intensity in Ne VIII (February 2, 1996)', yax='angstrom',xax='angstrom',charsize=1.5
write_png,'fullsun_NeVIII_02021996.png',tvrd(/true)
;stop

Free_lun,lun

;stop
fd_Neviii=rastor
;;;Save the 3D data cube as a sav file
;;save,fd_Neviii,filename='fd_Neviii.sav'


;for loop loop initializes and then runs through all fits files
;/silent doesnt print out
;checking for dimension of data and making sure it matches with header (Ex: CIV had larger fits files) ---> good to always check for dimension
;line 83 will print onto file (printf) if any of the data is too large/irrelevant
;/extract is for crvalues (want them in integer values) --->  FLOAT makes integer value
;don't want negative values ---> why sum doing +1024 (line 53)
;all reference pixels ---> we saw that CRVAL2 had 360/2 = 180.5 ---> have to put in correct position
;1 arcsecond in 1 pixel
;using array index (rastor) to make solar y ---> reference pixel changes (check the header to make sure what the reference pixel and what is the difference)

end

;[SAMPLE HEADER FILE FOR FEB 2, 1996]

;SIMPLE  =                    T / Written by IDL:  Tue Sep 11 15:24:13 2018
;BITPIX  =                  -32 /  IEEE single precision floating point
;NAXIS   =                    2 / Number of data axes
;NAXIS1  =                   25 /
;NAXIS2  =                  360 /
;EXTEND  =                    T / FITS data may contain extensions
;DATE    = '2018-09-07'         /  Creation date of FITS header
;FILENAME= 'sum_19960202_02051786_07701_12.fits' / Name of the FITS file
;SORIG   = 'GenuineIntel GNU/Linux' / Architecture and OS
;DATASRC = 'Final Distribution (CDROM)' / Data Source
;TELESCOP= 'SOHO    '           /
;INSTRUME= 'SUMER   '           /
;DATE_MOD= '2018-09-11T13:24:13.282' / Last modified
;ORIGIN  = 'SOHO MPS'           / Where Data is Produced
;OBS_SEQ = 'Raster  '           / Name of observing sequence
;DATE_OBS= '1996-02-02T02:05:17.864' / Beginning of Observation Date (UTC)
;DATE_END= '1996-02-02T02:05:25.365' / End of Observation Date (UTC)
;OBT_TIME=       1201917947.864 / Starting time of acquisition (TAI)
;OBT_END =       1201917955.365 / End time of acquisition (TAI)
;LEVEL   = '1       '           / Data processing level
;PRODLVL = 'L1      '           / Data calibration level
;DETECTOR= 'A       '           / Detector used
;EXPTIME =                 7.50 / Exposure time [s]
;IXWIDTH =              1.00000 / Image width ["]
;IYWIDTH =              300.000 / Image height ["]
;SLIT    = '<2> 1.0" * 300" centered' / SUMER slit
;INSTITUT= '        '           / Name of institutes for campaign
;CMP_TYPE= 'None    '           / Type of coordinated program
;CMP_NAME= 'None    '           / Name of campaign observation
;CMP_ID  =                    0 / Campaign number
;STUDY_ID=                  121 / Study number
;STUDY_NM= 'LEM POP06/2'        / Study name
;OBJECT  = 'FS      '           / Target
;SCI_OBJ = 'Full sun'           / The science objective
;SCI_SPEC= 'Full disk in Ne VIII 770' / The specified objective
;POPUDP  =                    6 / POP/UDP number
;OBS_PROG= 'CMP_SCHEME.SCL'     / Name of scientist program
;SCIENTIS= '< 19> Sys. Administrator' / Scientist responsible of POP/UDP
;FFONOFF = 'OFF     '           / On board flat field ON/OFF
;BINX    =                    1 / Binning X (1 - 1024)
;BINY    =                    1 / Binning Y (1 - 360)
;ROTCMP  =              0.00000 / Solar rotation compensation
;PIMGTYP = 'SINGLE  '           / Processing image Type (Single/Multiple)
;COMPRESS=                    5 / Compression method
;COMPAR1 =                   33 / Compression parameter 1
;COMPAR2 =                   20 / Compression parameter 2
;COMPAR3 =                    0 / Compression parameter 3
;DECOMP  =                    T / Decompression (F=No/T=Yes)
;ODEVCORR=                    T / Odd Even correction
;FLATCORR=                    T / Flatfield corrected on ground (F=No/T=Yes)
;FLATFILE= 'sumff_a_19960202_00422710_03.fits' / Used Flatfield Data
;FFMODDAT= '2016-10-26T06:32:56.212' / Creation/modification of Flatfield
;GEOMCORR=                    T / Geometry corrected on ground (F=No/T=Yes)
;DEADCORR=                    T / Dead time correction
;DCXDLEV =              9538.62 / Dead time corr XDL input value
;LGAINCOR=                    T / Local gain correction
;RADCORR =                    T / Radiometry calibration
;RADORDER=                    1 / Radiometry for Wavelength order
;AVARADO =                    1 / Available Radiometry orders
;WCSAXES =                    3 / Number of WCS axis
;CTYPE1  = 'WAVE    '           / Type of CRVAL1
;CUNIT1  = 'Angstrom'           / Axis unit along axis 1
;CRPIX1  =                   13 / Reference pixel along axis 1
;CRVAL1  =              770.159 / Value at reference pixel of axis 1
;CDELT1  =            0.0211100 / Axis increments along axis 1
;CTYPE2  = 'HPLT-TAN'           / Type of CRVAL2
;CUNIT2  = 'arcsec  '           / Axis unit along axis 2
;CRPIX2  =              180.500 / Reference pixel along axis 2
;CRVAL2  =             -947.749 / Value at reference pixel of axis 2
;CDELT2  =                    1 / Axis increments along axis 2
;CTYPE3  = 'HPLN-TAN'           / Type of CRVAL3
;CUNIT3  = 'arcsec  '           / Axis unit along axis 3
;CRPIX3  =                    1 / Reference pixel along axis 3
;CRVAL3  =             -623.779 / Value at reference pixel of axis 3
;CDELT3  =              1.00000 / Axis increments along axis 3 (Set to slit width
;IMGUNITS= 'W/sr/m^2/Angstroem' / Units for Image
;REFPIX  =                  830 / Original reference Pixel
;DETXSTRT=                  818 / Detector readout start / X
;DETXEND =                  842 / Detector readout end / X
;DETYSTRT=                    0 / Detector readout start / Y
;DETYEND =                  359 / Detector readout end / Y
;WAVEMIN =              769.906 / Minimum wavelength in image
;WAVEMAX =              770.412 / Maximum wavelength in image
;CORORB  =                    F / Orbitology corrected
;INS_X   =             -628.500 / Pointing of instrument / instrument X-axis
;INS_Y   =             -944.625 / Pointing of instrument / instrument Y-axis
;SOLAR_X =             -623.779 / Instrument pointing / solar X-axis
;SOLAR_Y =             -947.749 / Instrument pointing / solar Y-axis
;RASTYPE = 'RASTER  '           / Sequence type
;SC_X0   =              0.00000 / Spacecraft pointing / solar X-axis
;SC_Y0   =              0.00000 / Spacecraft pointing / solar Y-axis
;SC_ROLL =             0.285906 / Spacecraft roll angle relative to the solar coo
;SOLAR_P0=       -12.3397754705 / Solar angle P0 (degree)
;SOLAR_B0=       -6.07335233688 / Solar angle B0 (degree)
;GEIX_OBS=          3.76442e+08 / Geocentric equatorial inertial X
;GEIY_OBS=         -1.28532e+09 / Geocentric equatorial inertial Y
;GEIZ_OBS=         -5.79342e+08 / Geocentric equatorial inertial Z
;GSEX_OBS=          1.29288e+09 / Geocentric solar ecliptic X
;GSEY_OBS=         -6.76339e+08 / Geocentric solar ecliptic Y
;GSEZ_OBS=         -2.02509e+07 / Geocentric solar ecliptic Z
;GSMX_OBS=          1.29288e+09 / Geocentric solar magnetic X
;GSMY_OBS=         -6.08925e+08 / Geocentric solar magnetic Y
;GSMZ_OBS=         -2.95049e+08 / Geocentric solar magnetic Z
;HAEX_OBS=         -9.93361e+10 / Heliocentric Aries Ecliptic X
;HAEY_OBS=          1.07154e+11 / Heliocentric Aries Ecliptic Y
;HAEZ_OBS=         -1.57612e+07 / Heliocentric Aries Ecliptic Z
;DSUN_OBS=          1.46116e+11 / Distance from Sun
;HGLN_OBS=             0.229194 / Stonyhurst heliographic longitude
;HGLT_OBS=             -6.07335 / Stonyhurst heliographic latitude
;CRLN_OBS=              148.052 / Carrington heliographic longitude
;CRLT_OBS=             -6.07335 / Carrington heliographic latitude
;XCEN    =             -623.779 / Center of the instrument FOV / solar X-axis
;YCEN    =             -947.749 / Center of the instrument FOV / solar Y-axis
;ANGLE   =             0.285906 / Orientation of instrument FOV (degree)
;UDP_ID  =                  306 / Reference to observing program
;PROG_NM = 'UDP SYS CMP_SCH'    / Name of observing program
;T3TELE  =              52.1600 / Telescope (MC2) temperature (degree C)
;T3REAR  =              22.0000 / SUMER rear (MC3) temperature (degree C)
;T3FRONT =              23.6400 / SUMER front (MC4) temperature (degree C)
;T3SPACER=              22.5200 / SUMER spacer (MC6) temperature (degree C)
;MC2ENC  =                 2398 / SUMER MC2 Encoder Position
;MC3ENC  =                 2434 / SUMER MC3 Encoder Position
;MC4ENC  =                 1427 / SUMER MC4 Encoder Position
;MC6ENC  =                 3624 / SUMER MC6 Encoder Position
;MC8ENC  =               131072 / SUMER MC8 Encoder Position
;HKOTIME = '1996-02-02T02:05:20.342' / Time stamp of HK0 Record (UTC)
;
;        /SUMER Raw Image Header Entries
;
;SSYNC0  =                  235 / Header Sync Byte 0 (0xEB)
;SSYNC1  =                  144 / Header Sync Byte 1 (0x90)
;IMG_REC =                  128 / Header Record Id (0x80)
;SSTYPIMG=                   12 / SUMER image format index
;XSSTAI  =       1201917947.864 / Exposure start time (TAI)
;SSOPCNT =                   38 / Operations counter. 0 after switch on, increase
;SSPOPUDP=                    6 / POP/UDP number. 0 no POP/UDP executing
;SSIMGCNT=                    2 / Image counter, 0 at start of op +1 at spectrohe
;SSLOC   =                   19 / Scientist ID (after Jun 1996)
;SSTARGET=                   49 / Target as set by SYS_Operator
;SSFLDATE=                    0 / Never really used
;SSFLREQN=                    0 / Never really used
;SSREFPIX=                  830 / Reference pixel where lambda is on
;FREFPIX =                  830 / Ref pixel reduced by det B offset (2653)
;SSSTAT  =                   32 / Status
;SSEETRIG=                    0 / Explosive event trigger
;SSFF    =                    0 / Flatfield correction 1=on,0=off
;SSFF_T  = 'OFF     '           / Flatfield correction status as text
;SSDETTYP=                    1 / Current detector in use (1=A,2=B,3=RSC)
;SSDET_T = 'A       '           / Detector in use as text
;SSINTSTA=                    0 / IIF/TC spectrohelio interrupt status
;SSDETSTA=                  255 / Detector status
;SSSUNY  =               -10056 / SUMER coordinate Y [0.0625"]
;SSSUNZ  =                15114 / SUMER coordinate Z [0.0625"]
;P_SUNY  =             -628.500 / SUMER coordinate Y ["]
;P_SUNZ  =              944.625 / SUMER coordinate Z ["]
;SSEXPTIM=              7.50118 / Exposure time [s]
;SSIIDZ  =                    0 / Inter instrument sun z coordinate
;SSIIDY  =                    0 / Inter instrument sun y coordinate
;SSBPADDY=                    8 / Brightest pix addr Y (spec dir 0..1023)
;SSBPADDZ=                   51 / Brightest pix addr Z (spat dir 0..359)
;SSIMGTOT=                 4518 / Total counts in image
;SSROTCOM=              0.00000 / Rotation compensation time
;SSBPCNTS=                   20 / Brightest pixel counts
;SSACIMGC=                21961 / Accumulative image counter (0 after boot)
;SSSTEPN =                  663 / Number of raster steps (+ E->W, - W->E)
;SSSTEPSZ=                   -5 / Raster step size in motor steps (neg=Schmiersch
;SSSLITN =                    2 / Number of slit selected
;SSBINNY =                    1 / Binning factor spectral dimension
;SSBINNZ =                    1 / Binning factor spatial direction
;SSXCNT  =                 1365 / X event count (raw value)
;SSYCNT  =                 1365 / Y event count (raw value)
;X_EVCNT =              9538.62 / X event count (corrected for detector)
;Y_EVCNT =              9538.62 / Y event count (corrected for detector)
;S_MCPV  =                  218 / High voltage raw value
;S_MCPVF =              4835.00 / High voltage [V]
;SIMCPI  =                   96 / MCP current (Raw value)
;SIMCPIF =              29.4760 / MCP current [uA]
;SSMC2POS=                 4174 / MC2-Azimuth position [0.38"]
;SSMC3POS=                 3158 / MC3-Elev position [0.38"]
;SSMC4POS=                 -954 / MC4-Slit select position [half step]
;SSMC6POS=                18472 / MC6-Grat position [half step]
;SSMC8POS=                 9912 / MC8-Scan mirror position [half step]
;SSMCERR =                    0 / Motor controller error
;SSCOMPRM=                    5 / Compression method
;SSWAVEL =              770.159 / Wavelength at reference pixel (Angstroem)
;SSCOMPP1=                   33 / Compression parameter 1
;SSCOMPP2=                   20 / Compression parameter 2
;SSCOMPP3=                    0 / Compression parameter 3
;SSADMCNT=                  306 / Admin cnt (gives the database udp id)
;XSSMDU  =                    0 / Missing data in image 0=no,1=yes
;XSSQAC  =                    0 / Quality of image data 0=OK,1=NOTOK
;XSSFID  =                    0 / File id from TM file catalog
;XSSFPTR =              6732480 / Pointer to image position in bin file
;XSCDID  =                   85 / CD ID of TM file
;XSSEQID =                    6 / CD Sequence ID of TM file
;XSTMFILE=             60330202 / TM filename without ext
;PROCVERS=                  552 / SVN Version of write_sumerfits
;
;        / Data integrity check block
;
;CHECKSUM= 'WTanaTUlYTalaTUl'   / HDU checksum updated 2018-09-11T13:24:13
;DATASUM = '1719612052'         / data unit checksum updated 2018-09-11T13:24:13
;COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy
;COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H
;COMMENT The original raw image header is stored in the extention sumer-orig-raw-
;COMMENT The wavelength is stored low -> high
;COMMENT The data is stored with north up!
;COMMENT For documentation on SUMER FITS Data see:
;COMMENT Solar Physics June 2014, Volume 289, Issue 6, pp 2345–2376
;HISTORY Written with write_sumerfits Version:$Revision: 552 $
;HISTORY Applied decomp_method5 $Revision: 244 $
;HISTORY Applied update_fits $Revision: 561 $  7-Sep-2018 14:17:25.467
;HISTORY Applied DEADTIME_CORR Tue Sep 11 15:24:13 2018
;HISTORY Applied SUM_FLATFIELD with odd even corr Tue Sep 11 15:24:13 2018
;HISTORY Applied LOCAL_GAIN_CORR Tue Sep 11 15:24:13 2018
;HISTORY Applied SUM_FLATFIELD Tue Sep 11 15:24:13 2018
;HISTORY Applied DESTR_BILIN Tue Sep 11 15:24:13 2018
;HISTORY Applied RADIOMETRY Tue Sep 11 15:24:13 2018
;HISTORY Applied sumer_calib_fits Version:$Revision: 540 $
;HISTORY Calibration applied by SUMER
;END

