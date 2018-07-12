;-------------------------------------------------
;THIS IS THE ONLY ROUTINE YOU NEED TO BOTHER WITH.
;-------------------------------------------------
PRO satire_nlte
	;--------------------------------------------------------
	;CHANGE THIS TO THE DIRECTORY WHERE YOU HAVE THESE FILES.
	;--------------------------------------------------------
	path='/mnt/SSD/sim/satire/'
	;---------------------------------------------------
	;YOU PROVIDED TWO QS MODELS. I USED THE Q_KUR MODEL.
	;---------------------------------------------------
	qsn1_fn=path+'intensity_models/Q_kur'
	qsn2_fn=path+'intensity_models/Q_fal'
	fac1_fn=path+'intensity_models/F_fal'
	umb1_fn=path+'intensity_models/U_kur'
	pen1_fn=path+'intensity_models/P_kur'
	;-----------------------------------------------------------------------------------------------------------------------------------------
	;I SET THE RANGE OF BSAT TESTED AS 290G TO 300G. THE OPTIMAL VALUE IS 292G. AS YOU ADJUST THE SPECTRA YOU MIGHT HAVE TO CHANGE THIS RANGE.
	;THE ROUTINE WILL PRINT OUT THE OPTIMAL BSAT. IF YOU GET 290G OR 300G THEN YOU SHOULD EXTEND THIS RANGE.----------------------------------
	;-----------------------------------------------------------------------------------------------------------------------------------------

;	SATIRE,path,qsn1_fn,qsn2_fn,fac1_fn,umb1_fn,pen1_fn,290,300,/NLTE
;	SATIRE,path,qsn1_fn,qsn2_fn,fac1_fn,umb1_fn,pen1_fn,250,350,/NLTE
	SATIRE,path,qsn1_fn,qsn2_fn,fac1_fn,umb1_fn,pen1_fn,260,270,/NLTE
;	SATIRE,path,qsn1_fn,qsn2_fn,fac1_fn,umb1_fn,pen1_fn,266,266,/NLTE

	FIX_BSAT,path,free_p1,/NLTE

	SATIRE_SPEC,path,qsn1_fn,qsn2_fn,fac1_fn,umb1_fn,pen1_fn,free_p1,/NLTE

	RESTORE,path+'hmi_satire_date'
	RESTORE,path+'hmi_satire_tsi1'
	RESTORE,path+'hmi_satire_ssi1'

	tsi=INTERPOL(satire_tsi,satire_date,ROUND(satire_date))
	n=N_ELEMENTS(satire_date)	
	ssi=DBLARR(n,1221)	
	FOR i=0,1220 DO ssi[*,i]=INTERPOL(satire_ssi[*,i],satire_date,ROUND(satire_date))
	FOR i=0,n-1 DO ssi[i,*]=tsi[i]*ssi[i,*]/TSUM(wl,ssi[i,*])
	date=LONG(ROUND(satire_date))
	;---------------------------------------------------------------------------------------------
	;THIS IS THE FINAL OUTPUT FILE. THE RECONSTRUCTION BASED ON THE LTE SPECTRA IS SATIRE_LTE.SAV.
	;---------------------------------------------------------------------------------------------
	SAVE,date,wl,tsi,ssi,FILENAME=path+'satire_nlte.sav'

    print, 'done.'

END

PRO satire_lte

	path='/mnt/SSD/sim/satire/'

	qsn0_fn=path+'intensity_models/qs_1.dat'
	fac0_fn=path+'intensity_models/fac_5.dat'
	umb0_fn=path+'intensity_models/umbra_1.dat'
	pen0_fn=path+'intensity_models/penumbra_1.dat'

	SATIRE,path,qsn0_fn,fac0_fn,umb0_fn,pen0_fn,260,270,/LTE

	FIX_BSAT,path,free_p0,/LTE

	SATIRE_SPEC,path,qsn0_fn,fac0_fn,umb0_fn,pen0_fn,free_p0,/LTE

	RESTORE,path+'hmi_satire_date',/verbose
	RESTORE,path+'hmi_satire_tsi0',/verbose
	RESTORE,path+'hmi_satire_ssi0',/verbose
	tsi=INTERPOL(satire_tsi,satire_date,ROUND(satire_date))
	n=N_ELEMENTS(satire_date)		
	ssi=DBLARR(n,1221)	
	FOR i=0,1220 DO ssi[*,i]=INTERPOL(satire_ssi[*,i],satire_date,ROUND(satire_date))
	FOR i=0,n-1 DO ssi[i,*]=tsi[i]*ssi[i,*]/TSUM(wl,ssi[i,*])
	date=LONG(ROUND(satire_date))
	SAVE,date,wl,tsi,ssi,FILENAME=path+'satire_lte.sav'
END

PRO satire_spec,path,qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,free_p,lte=lte,nlte=nlte
	RESTORE,FILENAME=path+'hmi_satire_date'
	ff_fn=FILE_SEARCH(path+'fil_fac/*',COUNT=n)
	satire_ssi=DBLARR(n,1221)
	IF N_ELEMENTS(LTE)  NE 0 THEN RESTORE,FILENAME=path+'hmi_satire_tsi0'
	IF N_ELEMENTS(NLTE) NE 0 THEN RESTORE,FILENAME=path+'hmi_satire_tsi1'

	FOR i = 0, n - 1 DO begin

;        print, 'rel_irad_spec: date ', i + 1, ' out of ', n, ' dates'

        satire_ssi[i,*]=REL_IRAD_SPEC(qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,ff_fn[i],free_p,600.)

    endfor

	RDINTEN,qsn1_fn,d1,d2
	wl=d1[1,*]
	FOR i=0,n-1 DO satire_ssi[i,*]=satire_tsi[i]*satire_ssi[i,*]/TSUM(wl,satire_ssi[i,*])
	IF N_ELEMENTS(LTE)  NE 0 THEN SAVE,wl,satire_ssi,FILENAME=path+'hmi_satire_ssi0'
	IF N_ELEMENTS(NLTE) NE 0 THEN SAVE,wl,satire_ssi,FILENAME=path+'hmi_satire_ssi1'
END

PRO satire,path,qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,ll,ul,lte=lte,nlte=nlte
	IF N_ELEMENTS(LTE)  NE 0 THEN FILE_DELETE,path+'sat0/*'
	IF N_ELEMENTS(NLTE) NE 0 THEN FILE_DELETE,path+'sat1/*'
	ff_fn=FILE_SEARCH(path+'fil_fac/*',COUNT=n)
	satire_tsi=DBLARR(n)

	FOR free_p = ll, ul DO BEGIN

        print, 'free parameter: ', free_p, ' in range ', ll, ' ', ul

		FOR i = 0, n - 1 do begin

;            print, 'rel_irad: date ', i + 1, ' out of ', n, ' dates'

            satire_tsi[i]=REL_IRAD(qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,ff_fn[i],free_p,600.)

        endfor

		IF N_ELEMENTS(LTE)  NE 0 THEN SAVE,satire_tsi,FILENAME=path+'sat0/hmi_satire_tsi_'+STRCOMPRESS(free_p,/REMOVE_ALL)
		IF N_ELEMENTS(NLTE) NE 0 THEN SAVE,satire_tsi,FILENAME=path+'sat1/hmi_satire_tsi_'+STRCOMPRESS(free_p,/REMOVE_ALL)

	ENDFOR
END

PRO fix_bsat,path,free_p,lte=lte,nlte=nlte
	RESTORE,FILENAME=path+'hmi_satire_date'
	RESTORE,FILENAME=path+'measurements/tim.sav'
	IF N_ELEMENTS(lte)  NE 0 THEN list=FILE_SEARCH(path+'sat0/hmi_satire_tsi_*')
	IF N_ELEMENTS(nlte) NE 0 THEN list=FILE_SEARCH(path+'sat1/hmi_satire_tsi_*')
	;OVERLAP
	m=N_ELEMENTS(tim_date) & tim_cover=INTARR(m)
	n=N_ELEMENTS(satire_date) & satire_cover=INTARR(n)
	FOR i=0,m-1 DO BEGIN
		d=WHERE(ROUND(satire_date) EQ tim_date[i])
		IF d[0] NE -1 THEN BEGIN
			tim_cover[i]=1
			satire_cover[d]=1
		ENDIF
	ENDFOR
	tim_cover=WHERE(tim_cover EQ 1)
	tim=tim[tim_cover]
	satire_cover=WHERE(satire_cover EQ 1)
	;RMSD
	n=N_ELEMENTS(list)
	rmsd=DBLARR(n)
	conv=DBLARR(n)
	FOR i=0,n-1 DO BEGIN
		RESTORE,FILENAME=list[i]
		dsatire_tsi=INTERPOL(satire_tsi,satire_date,ROUND(satire_date))
		dsatire_tsi=dsatire_tsi[satire_cover]
		dsatire_tsi=1.+dsatire_tsi/1e6
		conv[i]=MEAN(tim)/MEAN(dsatire_tsi)
		dsatire_tsi=dsatire_tsi*conv[i]
		rmsd[i]=SQRT(MEAN((dsatire_tsi-tim)^2d))
	ENDFOR
	l=STRLEN(list[0])
	free_p=STRMID(list,l-3,3)
	d=MIN(rmsd,pos)
	free_p=free_p[pos]
	PRINT,'OPTIMAL BSAT: '+free_p
	IF N_ELEMENTS(LTE)  NE 0 THEN RESTORE,FILENAME=path+'sat0/hmi_satire_tsi_'+free_p
	IF N_ELEMENTS(NLTE) NE 0 THEN RESTORE,FILENAME=path+'sat1/hmi_satire_tsi_'+free_p
	satire_tsi=1.+satire_tsi/1e6
	satire_tsi=satire_tsi*conv[pos]
	IF N_ELEMENTS(LTE)  NE 0 THEN SAVE,satire_tsi,FILENAME=path+'hmi_satire_tsi0'
	IF N_ELEMENTS(NLTE) NE 0 THEN SAVE,satire_tsi,FILENAME=path+'hmi_satire_tsi1'
END

PRO calculate_julday
	ff_fn=FILE_SEARCH('/data/yeo/SATIRE-NLTE/fil_fac/*',COUNT=n)
	satire_date=DBLARR(n)
	l=STRLEN(ff_fn[0])
	FOR i=0,n-1 DO satire_date[i]=JULDAY(STRMID(ff_fn[i],l-26,2),STRMID(ff_fn[i],l-23,2),STRMID(ff_fn[i],l-31,4),$
                                         STRMID(ff_fn[i],l-20,2),STRMID(ff_fn[i],l-17,2),STRMID(ff_fn[i],l-14,2))
	SAVE,satire_date,FILENAME='/data/yeo/SATIRE-NLTE/hmi_satire_date'
END

PRO prepare_tim
	READCOL,'/data/yeo/SATIRE-NLTE/measurements/sorce_tsi_L3_c24h_latest.txt',d1,tim_date,d2,d3,tim,F='F,L,F,F,D',/SILENT
	d=WHERE(tim GT 0.)
	tim_date=tim_date[d]
	tim=tim[d]
	SAVE,tim_date,tim,FILENAME='/data/yeo/SATIRE-NLTE/measurements/tim.sav'
END

;----------------------------
;FACULAE FILLING FACTOR MODEL
;----------------------------
PRO ALPHA6,free_p,fac_hst,fac_alpha,sqrt=sqrt
	s=SIZE(fac_hst)
	fac_alpha=FLTARR(s[1],s[2])
	B=FINDGEN(s[2])*5.
	IF N_ELEMENTS(sqrt) EQ 0 THEN BEGIN
		FOR j=0,s[2]-1 DO IF B[j] GE free_p THEN fac_alpha[*,j]=fac_hst[*,j] ELSE fac_alpha[*,j]=fac_hst[*,j]*(B[j]/FLOAT(free_p))
	ENDIF ELSE BEGIN
		FOR j=0,s[2]-1 DO IF B[j] GE free_p THEN fac_alpha[*,j]=fac_hst[*,j] ELSE fac_alpha[*,j]=fac_hst[*,j]*SQRT(B[j]/FLOAT(free_p))
	ENDELSE
END

;---------------------------
;RELATIVE IRRADIANCE (TOTAL)
;---------------------------
FUNCTION REL_IRAD,qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,sav_fn,free_p,bcut,sqrt=sqrt
	;------------------------------------
	;LOAD INTENSITY MODELS & FILL FACTORS
	;------------------------------------
	RDINTEN,qsn1_fn,qsn1_int,qsn_mu
	RDINTEN,qsn2_fn,qsn2_int,qsn_mu
	RDINTEN,fac_fn,fac_int,fac_mu
	RDINTEN,umb_fn,umb_int,umb_mu
	RDINTEN,pen_fn,pen_int,pen_mu
	RESTORE,FILENAME=sav_fn
	;-------------------
	;RELATIVE IRRADIANCE
	;-------------------
	IF N_ELEMENTS(sqrt) EQ 0 THEN ALPHA6,free_p,fac_hst,fac_alpha ELSE ALPHA6,free_p,fac_hst,fac_alpha,/SQRT

	CALCINT,fac_alpha,umb_hst,pen_hst,mu_hst,qsn_mu,qsn1_int,qsn2_int,fac_int,umb_int,pen_int,image,image_qsn,fac_c,umb_c,pen_c,bcut

	ri=(TOTAL(image,/DOUBLE)/TOTAL(image_qsn,/DOUBLE)-1.d0)*1.d6
	RETURN,ri
END

;-------------------
;CALCULATE INTENSITY
;-------------------
PRO CALCINT,fac_ff,umb_ff,pen_ff,np,angles,qsn1_intlut,qsn2_intlut,fac_intlut,umb_intlut,pen_intlut,image,qsn_image,fac_c,umb_c,pen_c,bcut
;--------------
;INITIALIZATION
;--------------
	szm=SIZE(fac_ff)
	sza=SIZE(angles)

	image=DBLARR(szm[1])
	qsn_image=DBLARR(szm[1])

	qsn1_sumint=DBLARR(sza[1])
	qsn2_sumint=DBLARR(sza[1])

	fac_sumint=DBLARR(sza[1])
	umb_sumint=DBLARR(sza[1])
	pen_sumint=DBLARR(sza[1])

	c=DOUBLE(2.9979e8*1d9)					;SPEED OF LIGHT
	f=REFORM(c/(DOUBLE(qsn1_intlut[1,*])))	;FREQUENCY
;-------------------
;INTEGRATE INTENSITY
;-------------------
	FOR i=0,sza[1]-1 DO BEGIN
		d=REFORM(qsn1_intlut[i+3,*])*1d-3*f^2/c
		qsn1_sumint[i]=TSUM(qsn1_intlut[1,*],d)

		d=REFORM(qsn2_intlut[i+3,*])*1d-3*f^2/c
		qsn2_sumint[i]=TSUM(qsn2_intlut[1,*],d)

		d=REFORM(fac_intlut[i+3,*])*1d-3*f^2/c
		fac_sumint[i]=TSUM(fac_intlut[1,*],d)

		d=REFORM(umb_intlut[i+3,*])*1d-3*f^2/c
		umb_sumint[i]=TSUM(umb_intlut[1,*],d)

		d=REFORM(pen_intlut[i+3,*])*1d-3*f^2/c
		pen_sumint[i]=TSUM(pen_intlut[1,*],d)
	ENDFOR
;--------
;CONTRAST
;--------
	fac_c=fac_sumint/qsn2_sumint
	umb_c=umb_sumint/qsn1_sumint
	pen_c=pen_sumint/qsn1_sumint
;--------------
;INTERPOLATE MU
;--------------
	new_angles=FINDGEN(101)/100.

	new_qsn1_sumint=INTERPOL(qsn1_sumint,angles,new_angles)
	new_qsn2_sumint=INTERPOL(qsn2_sumint,angles,new_angles)

	new_fac_sumint=INTERPOL(fac_sumint,angles,new_angles)
	new_umb_sumint=INTERPOL(umb_sumint,angles,new_angles)
	new_pen_sumint=INTERPOL(pen_sumint,angles,new_angles)
;---------
;INTENSITY
;---------

	mu_max=.1

	cutoff = ROUND(mu_max/.01)

	IF bcut NE -1 THEN fac_ff[*,ROUND(bcut/5.):240]=0

	FOR i = cutoff, szm[1] - 1 DO BEGIN

        sf = new_qsn2_sumint[i] / new_qsn1_sumint[i]

		qsn_image[i] = np[i] * new_qsn2_sumint[i]

        image[i] = qsn_image[i] + total(fac_ff[i, *]) * (new_fac_sumint[i] - new_qsn2_sumint[i]) + $
                                        umb_ff[i]     * (new_umb_sumint[i] - new_qsn1_sumint[i]) * sf + $
                                        pen_ff[i]     * (new_pen_sumint[i] - new_qsn1_sumint[i]) * sf

    ENDFOR
END

;------------------------------
;RELATIVE IRRADIANCE (SPECTRAL)
;------------------------------
FUNCTION REL_IRAD_SPEC,qsn1_fn,qsn2_fn,fac_fn,umb_fn,pen_fn,sav_fn,free_p,bcut,sqrt=sqrt
	;------------------------------------
	;LOAD INTENSITY MODELS & FILL FACTORS
	;------------------------------------
	RDINTEN,qsn1_fn,qsn1_int,qsn_mu
	RDINTEN,qsn2_fn,qsn2_int,qsn_mu
	RDINTEN,fac_fn,fac_int,fac_mu
	RDINTEN,umb_fn,umb_int,umb_mu
	RDINTEN,pen_fn,pen_int,pen_mu
	RESTORE,FILENAME=sav_fn
	;-------------------
	;RELATIVE IRRADIANCE
	;-------------------
	IF N_ELEMENTS(sqrt) EQ 0 THEN ALPHA6,free_p,fac_hst,fac_alpha ELSE ALPHA6,free_p,fac_hst,fac_alpha,/SQRT
	CALCINT_SPEC,fac_alpha,umb_hst,pen_hst,mu_hst,qsn_mu,qsn1_int,qsn2_int,fac_int,umb_int,pen_int,image,image_qsn,fac_c,umb_c,pen_c,bcut
	ri=TOTAL(image,1,/DOUBLE)
	RETURN,ri
END

;-------------------
;CALCULATE INTENSITY
;-------------------
PRO CALCINT_SPEC,fac_ff,umb_ff,pen_ff,np,angles,qsn1_intlut,qsn2_intlut,fac_intlut,umb_intlut,pen_intlut,image,qsn_image,fac_c,umb_c,pen_c,bcut
;--------------
;INITIALIZATION
;--------------

    au    = 1.495985e+13 ; Astronomical unit (cm)
    R_sun = 6.9598e+10   ; solar radius (cm)

	szm = SIZE(fac_ff)      ; size new angles
	sza = SIZE(angles)      ; size old angles
	szl = SIZE(qsn1_intlut) ; size wavelengths

	image=DBLARR(szm[1],szl[2])
	qsn_image=DBLARR(szm[1],szl[2])

	qsn1_sumint=DBLARR(sza[1],szl[2])
	qsn2_sumint=DBLARR(sza[1],szl[2])

	fac_sumint=DBLARR(sza[1],szl[2])
	umb_sumint=DBLARR(sza[1],szl[2])
	pen_sumint=DBLARR(sza[1],szl[2])

    sf = DBLARR(szl[2])

	c=DOUBLE(2.9979e8*1d9)					;SPEED OF LIGHT
	f=REFORM(c/(DOUBLE(qsn1_intlut[1,*])))	;FREQUENCY
;-------------------
;INTEGRATE INTENSITY
;-------------------
	FOR i=0,sza[1]-1 DO BEGIN
		qsn1_sumint[i,*]=qsn1_intlut[i+3,*]*1d-3*f^2/c
		qsn2_sumint[i,*]=qsn2_intlut[i+3,*]*1d-3*f^2/c

		fac_sumint[i,*]=fac_intlut[i+3,*]*1d-3*f^2/c
		umb_sumint[i,*]=umb_intlut[i+3,*]*1d-3*f^2/c
		pen_sumint[i,*]=pen_intlut[i+3,*]*1d-3*f^2/c
	ENDFOR
;--------
;CONTRAST
;--------
	fac_c=fac_sumint/qsn2_sumint
	umb_c=umb_sumint/qsn1_sumint
	pen_c=pen_sumint/qsn1_sumint
;--------------
;INTERPOLATE MU
;--------------
	new_angles=FINDGEN(101)/100.

    ip = sqrt(1.0 - new_angles^2.0)

	new_qsn1_sumint=DBLARR(101,szl[2])
	new_qsn2_sumint=DBLARR(101,szl[2])

	new_fac_sumint=DBLARR(101,szl[2])
	new_umb_sumint=DBLARR(101,szl[2])
	new_pen_sumint=DBLARR(101,szl[2])

	FOR i=0,szl[2]-1 DO BEGIN
		new_qsn1_sumint[*,i]=INTERPOL(qsn1_sumint[*,i],angles,new_angles)
		new_qsn2_sumint[*,i]=INTERPOL(qsn2_sumint[*,i],angles,new_angles)

		new_fac_sumint[*,i]=INTERPOL(fac_sumint[*,i],angles,new_angles)
		new_umb_sumint[*,i]=INTERPOL(umb_sumint[*,i],angles,new_angles)
		new_pen_sumint[*,i]=INTERPOL(pen_sumint[*,i],angles,new_angles)
	ENDFOR
;---------
;INTENSITY
;---------
	mu_max = .1

	cutoff = ROUND(mu_max / .01)

	IF bcut NE -1 THEN fac_ff[*,ROUND(bcut/5.):240]=0

    sf[*] = 0.0d0

    for j = 0, szl[2] - 1 do begin

        qs2 = 0.0
        qs1 = 0.0

    	for i = cutoff, szm[1] - 2 do begin

            dOmega = !pi * (ip[i]^2.0 - ip[i + 1]^2.0) * (R_sun / au)^2.0

            if i ne szm[1] - 2 then begin

                I2 = (new_qsn2_sumint[i, j] + new_qsn2_sumint[i + 1, j]) / 2.0

                I1 = (new_qsn1_sumint[i, j] + new_qsn1_sumint[i + 1, j]) / 2.0

            endif else begin

                I2 = new_qsn2_sumint[szm[1] - 1, j]

                I1 = new_qsn1_sumint[szm[1] - 1, j]

            endelse

            qs2 = qs2 + I2 * dOmega

            qs1 = qs1 + I1 * dOmega

        endfor

        sf[j] = qs2 / qs1

        if sf[j] ne sf[j] then sf[j] = 1.0d0

    endfor

	FOR i = cutoff, szm[1] - 1 DO BEGIN

		qsn_image[i, *] = np[i] * new_qsn2_sumint[i, *]

        image[i, *] = qsn_image[i, *] + total(fac_ff[i, *]) * (new_fac_sumint[i, *] - new_qsn2_sumint[i, *]) + $
                                              umb_ff[i]     * (new_umb_sumint[i, *] - new_qsn1_sumint[i, *]) * sf[*] + $
                                              pen_ff[i]     * (new_pen_sumint[i, *] - new_qsn1_sumint[i, *]) * sf[*]

    ENDFOR

;print, 'szm[1] = ', szm[1]

;openw, lun1, 'test.txt', /get_lun

;for i = cutoff, szm[1] - 1 do begin

;    for j = 0, szl[2] - 1 do begin

;    printf, lun1, i, ' ', j, ' ', sf[i, j], ' ', image[i, j]

;    endfor

;endfor

;free_lun, lun1

END

PRO RDINTEN,filename,intens,angles
	angles=[1.,.9,.8,.7,.6,.5,.4,.3,.2,.1,.05]
	intens=DBLARR(14,1221)
	dummy=DBLARR(14)
	OPENR,1,filename
	FOR i=0,1220 DO BEGIN
		READF,1,dummy
		FOR j=0,13 DO intens(j,i)=dummy(j)
	ENDFOR
	CLOSE,1
END

FUNCTION TSUM,X,Y,IMIN,IMAX, NAN=NAN              ;Trapezoidal summation
;+
; NAME:
;       TSUM
; PURPOSE:
;       Trapezoidal summation of the area under a curve. 
; EXPLANATION:
;       Adapted from the procedure INTEG in the IUE procedure library.  
;
; CALLING SEQUENCE:
;       Result = TSUM(y)
;              or
;       Result = TSUM( x, y, [ imin, imax, /nan ] )  
; INPUTS:
;       x = array containing monotonic independent variable.  If omitted, then
;               x is assumed to contain the index of the y variable.
;               x = lindgen( N_elements(y) ).
;       y = array containing dependent variable y = f(x)
;
; OPTIONAL INPUTS:
;       imin = scalar index of x array at which to begin the integration
;               If omitted, then summation starts at x[0].
;       imax = scalar index of x value at which to end the integration 
;               If omitted then the integration ends at x[npts-1].
;       nan: If set cause the routine to check for occurrences of the IEEE 
;                 floating-point values NaN or Infinity in the input data. 
;                 Elements with the value NaN or Infinity are treated as missing data
;
; OUTPUTS:
;       result = area under the curve y=f(x) between x[imin] and x[imax].
;
; EXAMPLE:
;       IDL> x = [0.0,0.1,0.14,0.3] 
;       IDL> y = sin(x)
;       IDL> print,tsum(x,y)    ===>  0.0445843
;       
;       In this example, the exact curve can be computed analytically as 
;       1.0 - cos(0.3) = 0.0446635     
; PROCEDURE:
;       The area is determined of individual trapezoids defined by x[i],
;       x[i+1], y[i] and y[i+1].
;
;       If the data is known to be at all smooth, then a more accurate
;       integration can be found by interpolation prior to the trapezoidal
;       sums, for example, by the standard IDL User Library int_tabulated.pro.
; MODIFICATION HISTORY:
;       Written, W.B. Landsman, STI Corp. May 1986
;       Modified so X is not altered in a one parameter call Jan 1990
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Allow non-integer values of imin and imax  W. Landsman April 2001
;       Fix problem if only 1 parameter supplied W. Landsman June 2002
;       Added /nan keyword. Julio Castro/WL May 2014
;-
; Set default parameters
 On_error,2
 npar = N_params()
   
 if npar EQ 1 then begin
    npts = N_elements(x)
    yy = x
    xx = lindgen(npts)
    ilo = 0   & imin = ilo
    ihi = npts-1 & imax = ihi
 endif else begin

   if ( npar LT 3 ) then imin = 0
   npts = min( [N_elements(x), N_elements(y)] )
   if ( npar LT 4 ) then imax = npts-1
   ilo = long(imin)
   ihi = long(imax)
   xx = x[ilo:ihi]
   yy = y[ilo:ihi]
   npts = ihi - ilo + 1
 endelse   
; 
;  Remove NaN values
;
   if keyword_set(NaN) then begin 
   g = where(finite(yy),npts)
   yy = yy[g]
   xx = xx[g]
  endif          
;   
; Compute areas of trapezoids and sum result
;
  xdif = xx[1:*] - xx
  yavg =  ( yy[0:npts-2] + yy[1:npts-1] ) / 2.  
  sum = total( xdif*yavg ) 

; Now account for edge effects if IMIN or IMAX parameter are not integers

  hi = imax - ihi
  lo = imin - ilo
  if (ihi LT imax) then sum +=  (x[ihi+1]-x[ihi])*hi* $
              (y[ihi] + (hi/2.) *(y[ihi+1] - y[ihi]) )
  if (ilo LT imin) then sum -=  (x[ilo+1]-x[ilo])*lo* $
              (y[ilo] + (lo/2.) *(y[ilo+1] - y[ilo]) )
  return, sum

  end
