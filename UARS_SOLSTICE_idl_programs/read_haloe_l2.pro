; FILE: read_haloe_l2.pro
;
; SYNOPSIS: read_haloe_l2, 'filename', data, header, SWAP=swap
;
; DESCRIPTION:
;   This program will read in an entire HALOE level 2
;   data file, and save information to an array of data structures,
;   and a header structure.  Type `help, /str, header` and
;   `help, /str, data` to see contents of each structure, and
;   items within.
;
; AUTHOR: James Johnson (James.E.Johnson.1@gsfc.nasa.gov)
;
; HISTORY: Created Feb. 16, 2005

; --------

function l2_dat_rec

  return, { $

    head: { $
      label:'', nhead:0L, nhdlev:0L, nhdtyp:0L, $
      dates:0L, times:0L, datee:0L, timee:0L, mode:0L, nevent:0L, $
      sang:0.0, ainc:0.0, sz:0.0, zinc:0.0, $
      npts:0L, nrcrds:0L, iorb:0L, $
      salt:0.0, slat:0.0, slon:0.0, $
      nerror:lonarr(12), exosig:fltarr(12), sigval:fltarr(12), $
      erad90:0.0, erad30:0.0, erad6:0.0, $
      rdt:fltarr(4), stdev_rdt:fltarr(4), $
      filt_t:fltarr(4), stdev_filt_t:fltarr(4), $
      gc_t:fltarr(4), stdev_gc_t:fltarr(4), $
      beta:0.0, stlat:0.0, stlon:0.0, etlat:0.0, etlon:0.0, $
      evnlat:0.0, evnlon:0.0, evnvels:0.0, evnvela:0.0, method:intarr(4), $
      msisflag:0L, ch4_sat_z:0.0, ch4_sat_p:0.0, alt_gain:0.0, z_cirrus:0.0, $
      mch4:0L, evnstat:0L, ptflag:0L, smooth:intarr(12), indaero:intarr(12), $
      altlow:0.0, althigh:0.0, botexc:0.0, solextlo:0.0, apptoplo:0.0, $
      za_off_sun:0.0, ztrop:0.0, ptrop:0.0, ttrop:0.0, $
      idiflag:intarr(16) }, $
    data: replicate( { label:'', id:0L, npts:0L, value:fltarr(491) }, 281) $
  }

end

; --------

function l2_hdr_rec

  return, { $

    sfdu:'', $
    type:   { chfh:'', nwords:0L, $
              nhdlev2:0L, nhead:0L }, $
    comment:{ chfh:'', nwords:0L, $
              comment:'' }, $
    general:{ chfh:'', nwords:0L, $
              uars_day:0L, nl1events:0L, nret:0L, nskipped:0L }, $
    evnskip:{ chfh:'', nwords:0L, $
              skipped:lonarr(32) }, $
    sunset: { chfh:'', nwords:0L, $
              avglat_ss:0.0, avgvels_ss:0.0, avgvela_ss:0.0 }, $
    sunrise:{ chfh:'', nwords:0L, $
              avglat_sr:0.0, avgvels_sr:0.0, avgvela_sr:0.0 }, $
    evntype:{ chfh:'', nwords:0L, $
              evnt_type:strarr(32) }, $
    latinfo:{ chfh:'', nwords:0L, $
              sum_lat:fltarr(32) }, $
    loninfo:{ chfh:'', nwords:0L, $
              sum_lon:fltarr(32) }, $
    satvels:{ chfh:'', nwords:0L, $
              sum_vels:fltarr(32) }, $
    satvela:{ chfh:'', nwords:0L, $
              sum_vela:fltarr(32) }, $
    last:{ chfh:'', nwords:0L } }

end

; --------

PRO READ_F77_REC, unit, reclen, record, SWAP=swap

  reclen = 0L
  endlen = 0L

  READU, unit, reclen
  IF KEYWORD_SET(swap) THEN reclen = SWAP_ENDIAN(reclen)

  record = BYTARR(reclen)
  READU, unit, record
  READU, unit, endlen
  IF KEYWORD_SET(swap) THEN endlen = SWAP_ENDIAN(endlen)

  IF (reclen NE endlen) THEN reclen = -1

END

; --------

pro read_haloe_l2, drec, hdr, FILE=file, SWAP=swap

if not KEYWORD_SET(file) then begin

  file = ''
  print, format='("Please enter file name", $)
  read, file
  file = strtrim(file)

endif

  openr, unit, file, /get_lun

  hdr = l2_hdr_rec()

  ;Section 1: SFDU header

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  1 = SFDU Info
  hdr.sfdu = string(rec[0:reclen-1])

  ;Section 2: Summary section

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  2 = Header Type
  hdr.type.chfh         = string(rec[0:9])
  hdr.type.nwords       = long(rec,10,1)
  hdr.type.nhdlev2      = long(rec,14,1)
  hdr.type.nhead        = long(rec,18,1)
  IF KEYWORD_SET(swap) THEN hdr.type = SWAP_ENDIAN(hdr.type)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  3 = Comments
  hdr.comment.chfh      = string(rec[0:9])
  hdr.comment.nwords    = long(rec,10,1)
  hdr.comment.comment   = string(rec[14:reclen-1])
  IF KEYWORD_SET(swap) THEN hdr.comment = SWAP_ENDIAN(hdr.comment)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  4 = General Info
  hdr.general.chfh      = string(rec[0:9])
  hdr.general.nwords    = long(rec,10,1)
  hdr.general.uars_day  = long(rec,14,1)
  hdr.general.nl1events = long(rec,18,1)
  hdr.general.nret      = long(rec,22,1)
  hdr.general.nskipped  = long(rec,26,1)
  IF KEYWORD_SET(swap) THEN hdr.general = SWAP_ENDIAN(hdr.general)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  5 = Skipped Events
  hdr.evnskip.chfh      = string(rec[0:9])
  hdr.evnskip.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
    hdr.evnskip.nwords    = SWAP_ENDIAN(hdr.evnskip.nwords)
  hdr.evnskip.skipped  = long(rec,14,hdr.evnskip.nwords)
  IF KEYWORD_SET(swap) THEN $
    hdr.evnskip.skipped    = SWAP_ENDIAN(hdr.evnskip.skipped)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  6 = Avg Sunset Info
  hdr.sunset.chfh       = string(rec[0:9])
  hdr.sunset.nwords     = long(rec,10,1)
  hdr.sunset.avglat_ss  = float(rec,14,1)
  hdr.sunset.avgvels_ss = float(rec,18,1)
  hdr.sunset.avgvela_ss = float(rec,22,1)
  IF KEYWORD_SET(swap) THEN hdr.sunset = SWAP_ENDIAN(hdr.sunset)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  7 = Avg Sunrise Info
  hdr.sunrise.chfh      = string(rec[0:9])
  hdr.sunrise.nwords    = long(rec,10,1)
  hdr.sunrise.avglat_sr = float(rec,14,1)
  hdr.sunrise.avgvels_sr= float(rec,18,1)
  hdr.sunrise.avgvela_sr= float(rec,22,1)
  IF KEYWORD_SET(swap) THEN hdr.sunrise = SWAP_ENDIAN(hdr.sunrise)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  8 = Event Type Info
  hdr.evntype.chfh      = string(rec[0:9])
  hdr.evntype.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
  hdr.evntype.nwords = SWAP_ENDIAN(hdr.evntype.nwords)
  for i=0,hdr.evntype.nwords-1 do begin
    hdr.evntype.evnt_type[i] = string(rec[14+i*10:14+(i+1)*10-1])
  endfor

  read_f77_rec, unit, reclen, rec, swap=swap	; Record  9 = Latitude Info
  hdr.latinfo.chfh      = string(rec[0:9])
  hdr.latinfo.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
    hdr.latinfo.nwords = SWAP_ENDIAN(hdr.latinfo.nwords)
  hdr.latinfo.sum_lat[0:hdr.latinfo.nwords-1] = float(rec,14,hdr.latinfo.nwords)
  IF KEYWORD_SET(swap) THEN $
    hdr.latinfo.sum_lat = SWAP_ENDIAN(hdr.latinfo.sum_lat)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record 10 = Longitude Info
  hdr.loninfo.chfh      = string(rec[0:9])
  hdr.loninfo.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
    hdr.loninfo.nwords = SWAP_ENDIAN(hdr.loninfo.nwords)
  hdr.loninfo.sum_lon[0:hdr.latinfo.nwords-1] = float(rec,14,hdr.loninfo.nwords)
  IF KEYWORD_SET(swap) THEN $
    hdr.loninfo.sum_lon = SWAP_ENDIAN(hdr.loninfo.sum_lon)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record 11 = SC Vel to Sun Info
  hdr.satvels.chfh      = string(rec[0:9])
  hdr.satvels.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
    hdr.satvels.nwords = SWAP_ENDIAN(hdr.satvels.nwords)
  hdr.satvels.sum_vels[0:hdr.satvels.nwords-1]= float(rec,14,hdr.satvels.nwords)
  IF KEYWORD_SET(swap) THEN $
    hdr.satvels.sum_vels = SWAP_ENDIAN(hdr.satvels.sum_vels)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record 12 = SC Vel to Atm.Info
  hdr.satvela.chfh      = string(rec[0:9])
  hdr.satvela.nwords    = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN $
    hdr.satvela.nwords = SWAP_ENDIAN(hdr.satvela.nwords)
  hdr.satvela.sum_vela[0:hdr.satvela.nwords-1]= float(rec,14,hdr.satvela.nwords)
  IF KEYWORD_SET(swap) THEN $
    hdr.satvela.sum_vela = SWAP_ENDIAN(hdr.satvela.sum_vela)

  read_f77_rec, unit, reclen, rec, swap=swap	; Record 13 = End Header
  hdr.last.chfh         = string(rec[0:9])
  hdr.last.nwords       = long(rec,10,1)
  IF KEYWORD_SET(swap) THEN hdr.last = SWAP_ENDIAN(hdr.last)

  ;Loop through events

  drec = replicate( l2_dat_rec(), hdr.general.nret )

  n = 0

  while not eof(unit) do begin

    ;Event Header

    read_f77_rec, unit, reclen, rec, swap=swap

    drec[n].head.label        = string(rec[0:9])
    drec[n].head.nhead        = long(rec,10,1)
    drec[n].head.nhdlev       = long(rec,14,1)
    drec[n].head.nhdtyp       = long(rec,18,1)
    drec[n].head.dates        = long(rec,22,1)
    drec[n].head.times        = long(rec,26,1)
    drec[n].head.datee        = long(rec,30,1)
    drec[n].head.timee        = long(rec,34,1)
    drec[n].head.mode         = long(rec,38,1)
    drec[n].head.nevent       = long(rec,42,1)
    drec[n].head.sang         = float(rec,46,1)
    drec[n].head.ainc         = float(rec,50,1)
    drec[n].head.sz           = float(rec,54,1)
    drec[n].head.zinc         = float(rec,58,1)
    drec[n].head.npts         = long(rec,62,1)
    drec[n].head.nrcrds       = long(rec,66,1)
    drec[n].head.iorb         = long(rec,70,1)
    drec[n].head.salt         = float(rec,74,1)
    drec[n].head.slat         = float(rec,78,1)
    drec[n].head.slon         = float(rec,82,1)
    drec[n].head.nerror       = long(rec,86,12)
    drec[n].head.exosig       = float(rec,134,12)
    drec[n].head.sigval       = float(rec,182,12)
    drec[n].head.erad90       = float(rec,230,1)
    drec[n].head.erad30       = float(rec,234,1)
    drec[n].head.erad6        = float(rec,238,1)
    drec[n].head.rdt          = float(rec,242,4)
    drec[n].head.stdev_rdt    = float(rec,258,4)
    drec[n].head.filt_t       = float(rec,274,4)
    drec[n].head.stdev_filt_t = float(rec,290,4)
    drec[n].head.gc_t         = float(rec,306,4)
    drec[n].head.stdev_gc_t   = float(rec,322,4)
    drec[n].head.beta         = float(rec,338,1)
    drec[n].head.stlat        = float(rec,342,1)
    drec[n].head.stlon        = float(rec,346,1)
    drec[n].head.etlat        = float(rec,350,1)
    drec[n].head.etlon        = float(rec,354,1)
    drec[n].head.evnlat       = float(rec,358,1)
    drec[n].head.evnlon       = float(rec,362,1)
    drec[n].head.evnvels      = float(rec,366,1)
    drec[n].head.evnvela      = float(rec,370,1)
    drec[n].head.method       = fix(rec,374,4)
    drec[n].head.msisflag     = long(rec,382,1)
    drec[n].head.ch4_sat_z    = float(rec,386,1)
    drec[n].head.ch4_sat_p    = float(rec,390,1)
    drec[n].head.alt_gain     = float(rec,394,1)
    drec[n].head.z_cirrus     = float(rec,398,1)
    drec[n].head.mch4         = long(rec,402,1)
    drec[n].head.evnstat      = long(rec,406,1)
    drec[n].head.ptflag       = long(rec,410,1)
    drec[n].head.smooth       = fix(rec,414,12)
    drec[n].head.indaero      = fix(rec,438,12)
    drec[n].head.altlow       = float(rec,462,1)
    drec[n].head.althigh      = float(rec,466,1)
    drec[n].head.botexc       = float(rec,470,1)
    drec[n].head.solextlo     = float(rec,474,1)
    drec[n].head.apptoplo     = float(rec,478,1)
    drec[n].head.za_off_sun   = float(rec,482,1)
    drec[n].head.ztrop        = float(rec,486,1)
    drec[n].head.ptrop        = float(rec,490,1)
    drec[n].head.ttrop        = float(rec,494,1)
    drec[n].head.idiflag      = fix(rec,498,16)

    IF KEYWORD_SET(swap) THEN drec[n].head = SWAP_ENDIAN(drec[n].head)

  ;Event Data Records
    for i=0,drec[n].head.nrcrds-1 do begin

      read_f77_rec, unit, reclen, rec, swap=swap

      drec[n].data[i].label = string(rec[0:9])
      drec[n].data[i].id    = long(rec,10,1)
      IF KEYWORD_SET(swap) THEN $
        drec[n].data[i].id = SWAP_ENDIAN(drec[n].data[i].id)
      drec[n].data[i].npts  = long(rec,14,1)
      IF KEYWORD_SET(swap) THEN $
        drec[n].data[i].npts = SWAP_ENDIAN(drec[n].data[i].npts)

      if (drec[n].data[i].npts gt 0) then begin
        drec[n].data[i].value[0:drec[n].data[i].npts-1] = $
                                            float(rec,18,drec[n].data[i].npts)
      endif
      IF KEYWORD_SET(swap) THEN $
        drec[n].data[i].value = SWAP_ENDIAN(drec[n].data[i].value)

    endfor

    n = n+1

  endwhile

  skip: free_lun, unit

end

