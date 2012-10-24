function 

pro band_autocorr
;First, get the whole variability catalog of quasars on Stripe 82.

qsocat = qso_varcat()

;Next,index by thing_id, so we have separable catalogs of individual objects.

hist = histogram(qsocat.thing_id,reverse=ri,min=min(qsocat.thin_id))

;Count the number of distinct quasars.
nqso = total(hist gt 0)
nh = n_elements(hist)

for i = 0,nqso-1 do begin
    these_inds = ri[(ri[i]):(ri[i+1]-1)]
;   This contains all the observations for this object.
    these_qso = qsocat[these_inds]
;   Calculate the cross-correlation signal between all five bands as a
;   function of time.
    tvec = these_qso.mjd
    uflux = these_qso.psfflux[0]
    uivar = these_qso.psfflux_ivar[0]
    gflux = these_qso.psfflux[1]
    givar = these_qso.psfflux_ivar[1]
    rflux = these_qso.psfflux[2]
    rivar = these_qso.psfflux_ivar[2]
    iflux = these_qso.psfflux[3]
    iivar = these_qso.psfflux_ivar[3]
    zflux = these_qso.psfflux[4]
    zivar = these_Qso.psfflux_ivar[4]
    



endfor



end
