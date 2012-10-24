pro crosscorr_load,corr,lag,iqso,filter1,filter2
common Cmatrix_BLOCK,corr_matrix,lag_matrix,err_matrix
npts = n_elements(lag)
corr_matrix[0:(npts-1),iqso,filter1,filter2] = corr
lag_matrix[0:(npts-1),iqso,filter1,filter2] = lag
corr_matrix[*,iqso,filter2,filter1] = corr_matrix[*,iqso,filter1,filter2]
lag_matrix[ *,iqso,filter2,filter1] = lag_matrix[ *,iqso,filter1,filter2]

;err_matrix[*,iqso,filter1,filter2] = errvec
;err_matrix[*,iqso,filter2,filter1] = err_matrix[*,iqso,filter1,filter2]


end

pro band_crosscorr,nobj=nobj,qsocat,filename=filename
common Cmatrix_BLOCK,corr_matrix,lag_matrix
if ~keyword_set(filename) then filename = 'qso_preliminary_crosscorr'
;First, get the whole variability catalog of quasars on Stripe 82.

;qsocat = qso_varcat(nobj=nobj)
;restore,'qso_temp.sav'
;qsocat = mrdfits('../Data/qsocat_full.fits.gz',1)

qsocat = qsocat[where(qsocat.thing_id gt 0)]
;save,qsocat,file='qso_temp_full.sav'
print,'Successfully made QSO catalog.'
;Next,index by thing_id, so we have separable catalogs of individual objects.

hist = histogram(qsocat.thing_id,reverse=ri,min=min(qsocat.thing_id))

;Count the number of distinct quasars.
qso_use = where(hist gt 1,nqso)
nh = n_elements(hist)


;Make a grid to hold the final cross-correlations and lags, in each
;band. The lag grid is easy, since all the bands were taken
;simultaneously.
max_obs = max(hist)
lag_matrix = dblarr(max_obs,nqso,5,5)-1.
filterlist = ['u','g','r','i','z']
;For each possible filter pairing, we need a 200 by nqso matrix to
;hold the possible data.

corr_matrix = fltarr(max_obs,nqso,5,5)-1.
;A completely separate function is going to be needed to handle
;loading information into this beast.

err_matrix = fltarr(max_obs,nqso,5,5)-1.
print,'Making qso correlation matrix.'
for i = 0,nqso-1 do begin
    print,i,' quasars processed.'
    if ri[qso_use[i]] ne ri[qso_use[i]+1] then begin
        these_inds = ri[(ri[qso_use[i]]):(ri[qso_use[i]+1]-1)]
;   This contains all the observations for this object.
        if n_elements(these_inds) gt 1 then begin
            these_qso = qsocat[these_inds]
;   Calculate the cross-correlation signal between all five bands as a
;   function of time.
            tvec = these_qso.mjd
            tvec = tvec-min(tvec)
            for ifilt = 0,4 do begin
                for jfilt = 0,ifilt do begin
                    lagvec = findgen(n_elements(tvec))
                    cvec = c_correlate(these_qso.psfflux[ifilt],these_qso.psfflux[jfilt],lagvec)
                    crosscorr_load,cvec,tvec,i,ifilt,jfilt
                endfor
            endfor
        endif
    endif
endfor

;Now figure out how to organize the time gridding. We'll interpolate
;everything onto some final grid before averaging the cross-correlations.
min = 1.
max = 10000.
nbins = 8.
lagbins = min * 10.^(alog10( max / min ) * findgen(nbins+1) / (nbins))
;min = 500
;max = 2500
;lagbins = [500,1000,1500,2500,3000]

nbins = n_elements(lagbins)
lagbin_edges = [0,lagbins]
print,'All quasars processed. Proceeding to build binned correlation matrix.'
binned_corr_matrix = fltarr(nbins,5,5)
count_corr_matrix = fltarr(nbins,5,5)
for ifilt = 0,4 do begin
    for jfilt = 0,ifilt do begin
        for ibin = 0,nbins-2 do begin
;        Get all of the elements with lags in this range,
;        corresponding to this filter combination.
            these = where(lag_matrix[*,*,ifilt,jfilt] gt lagbin_edges[ibin] $
                          AND $
                          lag_matrix[*,*,ifilt,jfilt] le lagbin_edges[ibin+1],ct)

            if ct gt 0 then begin
                binned_corr_matrix[ibin,ifilt,jfilt] = total((corr_matrix[*,*,ifilt,jfilt])[these])/float(ct)
                binned_corr_matrix[ibin,jfilt,ifilt] = $
                  binned_corr_matrix[ibin,ifilt,jfilt]
                count_corr_matrix[ibin,ifilt,jfilt] = ct
                count_corr_matrix[ibin,jfilt,ifilt] = ct
            endif
        endfor
    endfor
endfor

psopen,filename,xsize=6,ysize=6,/inches
prepare_plots
for ifilt= 0,4 do begin
    for jfilt = 0,ifilt do begin
        plot,[0,lagbins],[1,binned_corr_matrix[*,ifilt,jfilt]],title=filterlist[ifilt]+' x '+filterlist[jfilt],xtitle='lag (days) ',ytitle='unweighted correlation'
    endfor
endfor
psclose
prepare_plots,/reset
stop
end
