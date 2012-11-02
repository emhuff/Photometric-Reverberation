pro estimator_test

;Generate a very long, very well-sampled lightcurve.

nobs = 50.
dt = 19.
tobs = findgen(nobs)*dt
;tobs =total(24.*randomu(seed,nobs),/cum)
tau=200.
tlag= 250.
psi_width = 80.
qso_lightcurve_sim,tobs,tlag=tlag,c=c,l=l,tau=tau,$
  transfer_sigma=psi_width,psi=psi

;This produces a visible peak in the cross-correlation.
;cl = c_correlate(c,l,findgen(nobs))
;ll = a_correlate(l,findgen(nobs))

;Add a small amount of white noise.
noise_amplitude = .050
noise = replicate(noise_amplitude^2,2*nobs)



;Now let us see if we can correctly back out the transfer function.
;We know that the lag is at t=100.
;Let's solve for psi in bins of width=20, ranging from t=0 to t=200.
w= 40.
tmax = 1000.
tmin = 0.
nbins = float(ceil((tmax-tmin)/w))
tpsi = [w/2+w*findgen(nbins)]


psi_in = exp(-(tpsi-tlag)^2/2./psi_width^2)
psi_in *= 0.
psi_in[0] = 1.0
;psi_in = interpol(psi,tobs[0:n_elements(psi)-1],tpsi)
;psi_in = psi_in/total(psi_in)

;Next, compute the covariance matrices and the optimal estimator and all.



niter= 20.
psi_avg = fltarr(nbins,niter)
for i=0L,niter-1 do begin
    undefine,c
    undefine,l
    qso_lightcurve_sim,tobs,tlag=tlag,c=c,l=l,tau=tau,$
      transfer_sigma=psi_width,psi=psi,sigma=1.0
    c += noise_amplitude*randomn(seed,nobs)
    l += noise_amplitude*randomn(seed,nobs)
    optimal_estimator,tobs,tpsi,[c,l],w,dt,sigma=1.0,mu=tau,psi_in=psi_in,$
      psi_out=psi_out_intermediate,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise
    optimal_estimator,tobs,tpsi,[c,l],w,dt,sigma=1.0,mu=tau,psi_in=psi_out_intermediate,$
      psi_out=psi_out,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise

    psi_avg[*,i] = -psi_out


    plot,tpsi,total(psi_avg,2)/float(niter)
    oplot,tobs,exp(-(tobs-tlag)^2/2./psi_width^2),color=200

endfor

prepare_plots,/color

filename = string(form='("../Plots/estimator_test.tlag.",I05,".ps")',long(tlag))

psopen,filename,/color
plot,tpsi,total(psi_avg,2)/10.,xtitle='t (days)',ytitle='!7 w !6'
oplot,tobs,exp(-(tobs-tlag)^2/2./psi_width^2),color=200
legend,['input','estimated'],line=0,color=[-1,200],/top,/right,box=0,charsize=2
psclose
prepare_plots,/reset



end

