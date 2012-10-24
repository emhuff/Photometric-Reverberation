function generate_qso_timeseries,lag,tlag=tlag,tau=tau,sigma=sigma
;Let's say for now that the only lines are in the i-band, and that the
;equivalent width is 2%.
qso_lightcurve_sim,lag,tlag=tlag,tau=tau,sigma=sigma,c=continuum,l=line
uflux = continuum
gflux = continuum
rflux = continuum
iflux = 0.98*continuum+0.02*line
zflux = continuum

;Add to each of these some white `measurement' noise with the same
;amplitude as the continuum signal.
amp = sigma
psfflux =transpose([[uflux],[gflux],[rflux],[iflux],[zflux]])
psfflux = psfflux + amp*randomn(seed,size(psfflux,/dim))

return,psfflux
end


pro qso_reverb_simulate,qsocat
COMMON qsocat_simulate_BLOCK,qsocat_orig,qsocat_sim
;if ~keyword_set(qsocat) then begin
;    qsocat_orig = mrdfits('../Data/qsocat_full.fits.gz',1)
;endif else 
qsocat_orig = qsocat
tlag= 1000.
qsocat_sim = qsocat_orig
nqso = n_elements(qsocat_orig)

;Make a unique list of thing_ids

thinglist = qsocat.thing_id
thinglist = thinglist[uniq(thinglist[sort(thinglist)])]
nthing = n_elements(thinglist)

for i = 0L,nthing-1 do begin
    these = where(qsocat.thing_id eq thinglist[i],ct)
    if ct gt 2 AND (max(qsocat[these].mjd) ne min(qsocat[these].mjd)) then begin
        lag = qsocat_sim[these].mjd
        qsocat_sim[these].psfflux = generate_qso_timeseries(lag,tlag=tlag)
    endif
endfor


band_crosscorr,qsocat_sim,filename='qso_crosscorr_sim.1000.ps'

stop

end
