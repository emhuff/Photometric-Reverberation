pro estimator_test

;Generate a very long, very well-sampled lightcurve.

nobs = 400.
dt = 7.3
tobs = findgen(nobs)*dt
;tobs =total(24.*randomu(seed,nobs),/cum)
tau=100.
tlag= 80.
psi_width = 20.
qso_lightcurve_sim,tobs,tlag=tlag,c=c,l=l,tau=tau,$
  transfer_sigma=psi_width,psi=psi

;This produces a visible peak in the cross-correlation.
;cl = c_correlate(c,l,findgen(nobs))
;ll = a_correlate(l,findgen(nobs))

;Add a small amount of white noise.
noise_amplitude = 1.0
noise = replicate(noise_amplitude,2*nobs)
c += noise_amplitude*randomn(seed,nobs)
l += noise_amplitude*randomn(seed,nobs)

;Now let us see if we can correctly back out the transfer function.
;We know that the lag is at t=100.
;Let's solve for psi in bins of width=20, ranging from t=0 to t=200.
w= 40.
tmax = 240.
tmin = 0.
nbins = float(ceil((tmax-tmin)/w))
tpsi = [w/2+w*findgen(nbins)]


psi_in = exp(-(tpsi-tlag)^2/2./psi_width^2)
psi_in *= 0.
psi_in[0] = 1.0
;psi_in = interpol(psi,tobs[0:n_elements(psi)-1],tpsi)
;psi_in = psi_in/total(psi_in)

;Next, compute the covariance matrices and the optimal estimator and all.

psi_avg = fltarr(nbins)

niter= 50.

for i=0L,niter-1 do begin
    undefine,c
    undefine,l
    qso_lightcurve_sim,tobs,tlag=tlag,c=c,l=l,tau=tau,$
      transfer_sigma=psi_width,psi=psi,sigma=1.0
    c += noise_amplitude*randomn(seed,nobs)
    l += noise_amplitude*randomn(seed,nobs)
    optimal_estimator,tobs,tpsi,[c,l],w,dt,sigma=1.0,mu=tau,psi_in=psi_in,$
      psi_out=psi_out,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise
    psi_avg += psi_out
    cl = c_correlate(c,l,findgen(nobs),/covariance)
    ll = a_correlate(l,findgen(nobs),/covariance)
    cc = a_correlate(c,findgen(nobs),/covariance)
    plot,tobs,LLCov[*,0],xr=[0,500]
    oplot,tobs,ll,color=1.5e7,thick=2
    

endfor
psi_avg = psi_avg/float(niter)

plot,tpsi,psi_avg


;optimal_estimator,tobs,tpsi,[c,l],w,sigma=1.0,mu=tau,psi_in=psi_out,$
;  psi_out=psi_final,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise


stop

end

