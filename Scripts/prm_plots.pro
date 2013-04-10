pro prm_plots


;Make well-sampled quasar lightcurves.
restore,'example_mjd.sav'
dt = 5.
tobs = findgen(200)*dt
;mjd = mjd-min(mjd)
;tobs = mjd

mu = 50.
tlag = 45.0
sigma_psi = 10.


nobs = n_elements(tobs)
tpsi = findgen(50)*4.
psi = exp(-(tpsi-tlag)^2/2./sigma_psi^2)/sqrt(2*!Pi*sigma_psi^2)
psi = 100.*psi/int_tabulated(tpsi,psi)
noise = replicate(0.1,2*nobs)
nqso = 100.
mixing = 0.9
psi_all = fltarr(n_elements(psi),nqso)
photo_psi = psi_all
variance = psi_all

for i = 0,nqso-1 do begin
   print,i
   data = drw_sim_covar(tobs,tpsi,psi,/reset,cond=cond,sigma=1.,mu=mu,noise=noise)
   continuum = data[0:(nobs-1)]
   line = data[nobs:(2*nobs-1)]
   data_photo = [continuum,0.9*continuum+0.1*line]
;   psopen,'example_lightcurve_well-sampled.ps',xsize=6,ysize=6,/inches,/enc,/color
;   prepare_plots,/color
;   plot,tobs,continuum+3,xtitle='t (days)',ytitle='flux (arbitrary)',yr=[-8,8]
;   oplot,tobs,line/10.-3,color=150
;   oplot,[tobs[100],tobs[100]],[-10,continuum[100]+3],color=200
;   oplot,[tobs[100],tobs[100]]+tlag,[-10,interpol(line/10.-3,tobs,tobs[100]+tlag)],color=200
;   prepare_plots,/reset
;   psclose
   
   w = tobs[1]-tobs[0]
   
   optimal_estimator,tobs,tpsi,data,w,1.0,sigma=1.0,mu=tau,psi_in=psi,psi_out=psi_out,fisher=fisher,noise=noise
   optimal_estimator,tobs,tpsi,data_photo,w,1.0,sigma=1.0,mu=tau,psi_in=psi,psi_out=psi_out_photo,noise=noise
   psi_all[*,i] = psi_out
   photo_psi[*,i] = psi_out_photo
   variance[*,i] = diag_matrix(inverT(fisher,/double))
endfor
ploterror,tpsi,total(psi_all/variance,2)/total(variance,2),sqrt(total(variance,2)/100.)
stop

end
