; NAME:
;  fsc_NeVIII_06161996.PRO
;
; PURPOSE:
;  To construct full sun 3D data cube from the slit scan of NE VIII

; CALLING SEQUENCE:
;  FULLSUN_NeVIII_06161996

; INPUTS:
;  NONE

; OPTIONAL INPUTS:
;  NONE
;
; OUTPUTS:
;  *sav file with the 3D data cube (also saves to drive)
;  *png. file of Full Disk Peak Intensity of NeVIII [June 16, 1996] (also saves to drive)

; MODIFICATION HISTORY:
;  Written by: Suman Panda 14 Jun 2021 
;  Last Modified: Alyssa Johnson 27 Jun 2021 
; 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro fullsun_NeVIII_06161996

OpenW,lun,'fulldisk_NeVIII_06161996.txt',/Get_Lun	

rastor=fltarr(1,3000,1000)
rastor_step=0
ys0=0
ys=0
xs=0
xss=0
strip_no=1

drive='/disk/data/ajohnson/CIV_NeVIII_fitsfiles/1996/06/16/'
fits_files=file_search(drive,'*.fits',count=count)	
	
fps=3
ovid=idlffvideowrite('fullsun_NEVIII_06161996.mp4')
window,xs=800,ys=800
image=tvrd(/true)
s=size(image)
xsize=s[2]
ysize=s[3]
vidstream=ovid.addvideostream(xsize,ysize,fps)
 
 
for i= 0, count-1 do begin
	data= readfits(fits_files(i),Header,/silent)
	
	dimention=size(data)
	if dimention(1) eq 1 and dimention(2) eq 360 then begin
		xs0=xss
		CRVAL2_str=Header(59)
		substr1=strsplit(crval2_str,' ',/extract)
		CRVAL2=FLOAT(substr1(2))

		CRVAL3_str=Header(64)
		substr2=strsplit(crval3_str,' ',/extract)
		CRVAL3=FLOAT(substr2(2))
			
		xss=ROUND((CRVAL3+1024))/3
		yss=ROUND(CRVAL2+1500-60.5)
		
		if abs(ys0-yss) ge 250 then begin
			ys=yss
			xs=xss;-abs(xs0-xs)		
			strip_no=strip_no+1
			window,xs=800,ys=800
			plot_image,sqrt(max(rastor,dimension=1)),title='Full Disk Peak Intensity in Ne VIII (June 16, 1996)', xtitle='arcseconds',ytitle='arcseconds'
			image2=tvrd(/true) ;grabs window currently displayed, if you put before window nothing would be grabbed
      time=ovid.put(vidstream,image2) ;error was from vidstream and image NOT being the same size
      ;cursor,x,y	
		endif				
		 	
		rastor[0,yss:yss+299,xss]=data[0,60:359]
		ys0=yss
		
		if xs ne xss then rastor[0,yss:yss+299,xs]=data[0,60:359]

		if (strip_no MOD 2) eq 1 then xs=xss-1 else xs=xss+1
	endif else begin
		;printf,lun,fits_files(i)
	endelse
	
endfor
	window,xs=800,ys=800
	tvscl,sqrt(max(rastor,dimension=1))
	cursor,x,y

peak_intensity = max(rastor,dimension=1)
peak_intensity = TRANSPOSE(peak_intensity)
window,xs=800,ys=800

plot_image,sqrt(max(rastor,dimension=1)),title='Full Disk Peak Intensity in Ne VIII (June 16, 1996)', xtitle='arcseconds',ytitle='arcseconds',charsize=1.5
write_png,'fullsun_NeVIII_06161996.png',tvrd(/true)

Free_lun,lun

ovid.cleanup 

end


;;;;;;;;;;;;;;;;;
;CODE NOTES
;;;;;;;;;;;;;;;;;
;this keeps track of any files that have dimention other than the 25x360
;OpenW,lun,'fulldisk_NeVIII_06071996.txt',/Get_Lun
;this initializes the array
;the first entry is for NAXIS2
;the second and third entries you get from using CALLING SEQUENCE plot_xy_position_06071996 aka count value
;X and y value of where they are on the Sun
;for i= 0, count-1 do begin ;adding,3 goes in steps of three
;  rastor[*,ys:ys+299,xs]=data[0,50:359]
;print,'i=',i
;fits_files=file_search(drive,'*.fits',count=count)
;print,count
;stop	
;peak_intensity = max(rastor,dimension=1)
;peak_intensity = TRANSPOSE(peak_intensity)
;window,xs=1000,ys=1000
;plot_image,sqrt(congrid(peak_intensity,1000,1000))
;Free_lun,lun
;stop
;fd_Neviii=rastor
;;;Save the 3D data cube as a sav file
;;save,fd_Neviii,filename='fd_Neviii.sav'
;end


;;;;;;;;;;;;;;;;;;
;SUMAN NOTES:
;;;;;;;;;;;;;;;;;;
;(1) for loop initializes and then runs through all fits files
;(2) /silent doesnt print out
;(3) checking for dimension of data and making sure it matches with header (Ex: CIV had larger fits files) ---> good to always check for dimension
;(4) line 83 will print onto file (printf) if any of the data is too large/irrelevant
;(5)  /extract is for crvalues (want them in integer values) --->  FLOAT makes integer value
;(6)  don't want negative values ---> why sum doing +1024 (line 53)
;(7)  all reference pixels ---> we saw that CRVAL2 had 360/2 = 180.5 ---> have to put in correct position
;(8)  1 arcsecond in 1 pixel
;(9) using array index (rastor) to make solar y ---> reference pixel changes (check the header to make sure what the reference pixel and what is the difference)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;[SAMPLE HEADER FILE FOR ----> June 16, 1996]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SIMPLE  =                    T / Written by IDL:  Sat Sep  8 05:53:16 2018
;BITPIX  =                    8 / Number of bits per data pixel
;NAXIS   =                    2 / Number of data axes
;NAXIS1  =                    1 /
;NAXIS2  =                  360 /
;EXTEND  =                    T / FITS data may contain extensions
;DATE    = '2018-09-08'         /  Creation date of FITS header
;FILENAME= 'sum_19960616_11313010_07877_20.fits' / Name of the FITS file
;SORIG   = 'GenuineIntel GNU/Linux' / Architecture and OS
;DATASRC = 'Final Distribution (CDROM)' / Data Source
;TELESCOP= 'SOHO    '           /
;INSTRUME= 'SUMER   '           /
;DATE_MOD= '2018-09-08T03:53:16.036' / Last modified
;ORIGIN  = 'SOHO MPS'           / Where Data is Produced
;OBS_SEQ = 'Raster  '           / Name of observing sequence
;DATE_OBS= '1996-06-16T11:31:30.101' / Beginning of Observation Date (UTC)
;DATE_END= '1996-06-16T11:31:32.103' / End of Observation Date (UTC)
;OBT_TIME=       1213615920.101 / Starting time of acquisition (TAI)
;OBT_END =       1213615922.103 / End time of acquisition (TAI)
;LEVEL   = '0       '           / Data processing level
;PRODLVL = 'LZ      '           / Data calibration level
;DETECTOR= 'A       '           / Detector used
;EXPTIME =                 2.00 / Exposure time [s]
;IXWIDTH =              1.00000 / Image width ["]
;IYWIDTH =              300.000 / Image height ["]
;SLIT    = '<2> 1.0" * 300" centered' / SUMER slit
;INSTITUT= '        '           / Name of institutes for campaign
;CMP_TYPE= 'None    '           / Type of coordinated program
;CMP_NAME= 'None    '           / Name of campaign observation
;CMP_ID  =                    0 / Campaign number
;STUDY_ID=                  344 / Study number
;STUDY_NM= 'HAS FS 1548/770'    / Study name
;OBJECT  = 'FS      '           / Target
;SCI_OBJ = 'Full sun'           / The science objective
;SCI_SPEC= 'Full Sun in C IV, O IV, Ne VIII HAS' / The specified objective
;POPUDP  =                    6 / POP/UDP number
;OBS_PROG= 'POP_K.SCL'          / Name of scientist program
;SCIENTIS= '<  9> Philippe Lemaire' / Scientist responsible of POP/UDP
;FFONOFF = 'OFF     '           / On board flat field ON/OFF
;BINX    =                    1 / Binning X (1 - 1024)
;BINY    =                    3 / Binning Y (1 - 360)
;ROTCMP  =              0.00000 / Solar rotation compensation
;PIMGTYP = 'SINGLE  '           / Processing image Type (Single/Multiple)
;COMPRESS=                   16 / Compression method
;COMPAR1 =                    0 / Compression parameter 1
;COMPAR2 =                    0 / Compression parameter 2
;COMPAR3 =                    0 / Compression parameter 3
;DECOMP  =                    T / Decompression (F=No/T=Yes)
;FLATCORR=                    F / Flatfield corrected on ground (F=No/T=Yes)
;GEOMCORR=                    F / Geometry corrected on ground (F=No/T=Yes)
;WCSAXES =                    3 / Number of WCS axis
;CTYPE1  = 'WAVE    '           / Type of CRVAL1
;CUNIT1  = 'Angstrom'           / Axis unit along axis 1
;CRPIX1  =                    1 / Reference pixel along axis 1
;CRVAL1  =              787.700 / Value at reference pixel of axis 1
;CDELT1  =            0.0210955 / Axis increments along axis 1
;CTYPE2  = 'HPLT-TAN'           / Type of CRVAL2
;CUNIT2  = 'arcsec  '           / Axis unit along axis 2
;CRPIX2  =              60.5000 / Reference pixel along axis 2
;CRVAL2  =             -946.391 / Value at reference pixel of axis 2
;CDELT2  =                    3 / Axis increments along axis 2
;CTYPE3  = 'HPLN-TAN'           / Type of CRVAL3
;CUNIT3  = 'arcsec  '           / Axis unit along axis 3
;CRPIX3  =                    1 / Reference pixel along axis 3
;CRVAL3  =             -622.579 / Value at reference pixel of axis 3
;CDELT3  =              1.00000 / Axis increments along axis 3 (Set to slit width
;REFPIX  =                  121 / Original reference Pixel
;DETXSTRT=                  121 / Detector readout start / X
;DETXEND =                  121 / Detector readout end / X
;WAVEMIN =              787.700 / Minimum wavelength in image
;WAVEMAX =              787.700 / Maximum wavelength in image
;CORORB  =                    F / Orbitology corrected
;INS_X   =             -624.688 / Pointing of instrument / instrument X-axis
;INS_Y   =             -945.000 / Pointing of instrument / instrument Y-axis
;SOLAR_X =             -622.579 / Instrument pointing / solar X-axis
;SOLAR_Y =             -946.391 / Instrument pointing / solar Y-axis
;RASTYPE = 'RASTER  '           / Sequence type
;SC_X0   =         -3.16734e-05 / Spacecraft pointing / solar X-axis
;SC_Y0   =            0.0142033 / Spacecraft pointing / solar Y-axis
;SC_ROLL =             0.127770 / Spacecraft roll angle relative to the solar coo
;SOLAR_P0=       -9.05414149877 / Solar angle P0 (degree)
;SOLAR_B0=        1.20321130753 / Solar angle B0 (degree)
;GEIX_OBS=          2.23315e+08 / Geocentric equatorial inertial X
;GEIY_OBS=          1.59730e+09 / Geocentric equatorial inertial Y
;GEIZ_OBS=          5.42632e+08 / Geocentric equatorial inertial Z
;GSEX_OBS=          1.69344e+09 / Geocentric solar ecliptic X
;GSEY_OBS=         -9.49176e+07 / Geocentric solar ecliptic Y
;GSEZ_OBS=         -1.37523e+08 / Geocentric solar ecliptic Z
;GSMX_OBS=          1.69344e+09 / Geocentric solar magnetic X
;GSMY_OBS=         -1.16572e+08 / Geocentric solar magnetic Y
;GSMZ_OBS=         -1.19720e+08 / Geocentric solar magnetic Z
;HAEX_OBS=         -1.13252e+10 / Heliocentric Aries Ecliptic X
;HAEY_OBS=         -1.49870e+11 / Heliocentric Aries Ecliptic Y
;HAEZ_OBS=         -1.44553e+08 / Heliocentric Aries Ecliptic Z
;DSUN_OBS=          1.50297e+11 / Distance from Sun
;HGLN_OBS=              0.00000 / Stonyhurst heliographic longitude
;HGLT_OBS=              1.20321 / Stonyhurst heliographic latitude
;CRLN_OBS=              161.001 / Carrington heliographic longitude
;CRLT_OBS=              1.20321 / Carrington heliographic latitude
;XCEN    =             -622.579 / Center of the instrument FOV / solar X-axis
;YCEN    =             -946.391 / Center of the instrument FOV / solar Y-axis
;ANGLE   =             0.127770 / Orientation of instrument FOV (degree)
;UDP_ID  =                  720 / Reference to observing program
;PROG_NM = 'POP 06  '           / Name of observing program
;T3TELE  =              51.2000 / Telescope (MC2) temperature (degree C)
;T3REAR  =              19.9200 / SUMER rear (MC3) temperature (degree C)
;T3FRONT =              20.0000 / SUMER front (MC4) temperature (degree C)
;T3SPACER=              23.0800 / SUMER spacer (MC6) temperature (degree C)
;MC2ENC  =                 2400 / SUMER MC2 Encoder Position
;MC3ENC  =                 2074 / SUMER MC3 Encoder Position
;MC4ENC  =                 1434 / SUMER MC4 Encoder Position
;MC6ENC  =                 3382 / SUMER MC6 Encoder Position
;MC8ENC  =               131072 / SUMER MC8 Encoder Position
;HKOTIME = '1996-06-16T11:31:24.019' / Time stamp of HK0 Record (UTC)
;
;        /SUMER Raw Image Header Entries
;
;SSYNC0  =                  235 / Header Sync Byte 0 (0xEB)
;SSYNC1  =                  144 / Header Sync Byte 1 (0x90)
;IMG_REC =                  128 / Header Record Id (0x80)
;SSTYPIMG=                   20 / SUMER image format index
;XSSTAI  =       1213615920.101 / Exposure start time (TAI)
;SSOPCNT =                  232 / Operations counter. 0 after switch on, increase
;SSPOPUDP=                    6 / POP/UDP number. 0 no POP/UDP executing
;SSIMGCNT=                    3 / Image counter, 0 at start of op +1 at spectrohe
;SSLOC   =                   22 / Scientist ID (after Jun 1996)
;SSTARGET=                    0 / Target as set by SYS_Operator
;SSFLDATE=                    0 / Never really used
;SSFLREQN=                    0 / Never really used
;SSREFPIX=                  121 / Reference pixel where lambda is on
;FREFPIX =                  121 / Ref pixel reduced by det B offset (2653)
;SSSTAT  =                   32 / Status
;SSEETRIG=                    0 / Explosive event trigger
;SSFF    =                    0 / Flatfield correction 1=on,0=off
;SSFF_T  = 'OFF     '           / Flatfield correction status as text
;SSDETTYP=                    1 / Current detector in use (1=A,2=B,3=RSC)
;SSDET_T = 'A       '           / Detector in use as text
;SSINTSTA=                    0 / IIF/TC spectrohelio interrupt status
;SSDETSTA=                  255 / Detector status
;SSSUNY  =                -9995 / SUMER coordinate Y [0.0625"]
;SSSUNZ  =                15120 / SUMER coordinate Z [0.0625"]
;P_SUNY  =             -624.688 / SUMER coordinate Y ["]
;P_SUNZ  =              945.000 / SUMER coordinate Z ["]
;SSEXPTIM=              2.00218 / Exposure time [s]
;SSIIDZ  =                    0 / Inter instrument sun z coordinate
;SSIIDY  =                    0 / Inter instrument sun y coordinate
;SSBPADDY=                   35 / Brightest pix addr Y (spec dir 0..1023)
;SSBPADDZ=                   10 / Brightest pix addr Z (spat dir 0..359)
;SSIMGTOT=                    9 / Total counts in image
;SSROTCOM=              0.00000 / Rotation compensation time
;SSBPCNTS=                    1 / Brightest pixel counts
;SSACIMGC=                 1862 / Accumulative image counter (0 after boot)
;SSSTEPN =                  414 / Number of raster steps (+ E->W, - W->E)
;SSSTEPSZ=                   -8 / Raster step size in motor steps (neg=Schmiersch
;SSSLITN =                    2 / Number of slit selected
;SSBINNY =                    1 / Binning factor spectral dimension
;SSBINNZ =                    3 / Binning factor spatial direction
;SSXCNT  =                  296 / X event count (raw value)
;SSYCNT  =                  298 / Y event count (raw value)
;X_EVCNT =              2068.45 / X event count (corrected for detector)
;Y_EVCNT =              2082.42 / Y event count (corrected for detector)
;S_MCPV  =                  238 / High voltage raw value
;S_MCPVF =              5279.00 / High voltage [V]
;SIMCPI  =                  116 / MCP current (Raw value)
;SIMCPIF =              35.7160 / MCP current [uA]
;SSMC2POS=                 4170 / MC2-Azimuth position [0.38"]
;SSMC3POS=                 2540 / MC3-Elev position [0.38"]
;SSMC4POS=                 -954 / MC4-Slit select position [half step]
;SSMC6POS=                18664 / MC6-Grat position [half step]
;SSMC8POS=                10054 / MC8-Scan mirror position [half step]
;SSMCERR =                    0 / Motor controller error
;SSCOMPRM=                   16 / Compression method
;SSWAVEL =              787.700 / Wavelength at reference pixel (Angstroem)
;SSCOMPP1=                    0 / Compression parameter 1
;SSCOMPP2=                    0 / Compression parameter 2
;SSCOMPP3=                    0 / Compression parameter 3
;SSADMCNT=                  720 / Admin cnt (gives the database udp id)
;XSSMDU  =                    0 / Missing data in image 0=no,1=yes
;XSSQAC  =                    0 / Quality of image data 0=OK,1=NOTOK
;XSSFID  =                    0 / File id from TM file catalog
;XSSFPTR =             37541640 / Pointer to image position in bin file
;XSCDID  =                   63 / CD ID of TM file
;XSSEQID =                    3 / CD Sequence ID of TM file
;XSTMFILE=             61680101 / TM filename without ext
;PROCVERS=                  552 / SVN Version of write_sumerfits
;
;        / Data integrity check block
;
;CHECKSUM= 'BXRXEWQWBWQWBWQW'   / HDU checksum updated 2018-09-08T03:53:16
;DATASUM = '0       '           / data unit checksum updated 2018-09-08T03:53:16
;COMMENT FITS (Flexible Image Transport System) format is defined in 'Astronomy
;COMMENT and Astrophysics', volume 376, page 359; bibcode 2001A&A...376..359H
;COMMENT The original raw image header is stored in the extention sumer-orig-raw-
;COMMENT The wavelength is stored low -> high
;COMMENT The data is stored with north up!
;COMMENT For documentation on SUMER FITS Data see:
;COMMENT Solar Physics June 2014, Volume 289, Issue 6, pp 2345–2376
;HISTORY Written with write_sumerfits Version:$Revision: 552 $
;HISTORY Applied decomp_method16 $Revision: 244 $
;HISTORY Applied update_fits $Revision: 561 $  8-Sep-2018 03:53:16.036
;END
