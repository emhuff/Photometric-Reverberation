pro estimator_test,tlag=tlag,nobs=nobs,w=w

;Generate a very long, very well-sampled lightcurve.

;if ~keyword_set(nobs) then nobs = 50.
;dt = 1000/float(nobs)
;tobs = findgen(nobs)*dt
restore,'example_mjd.sav'
tobs = [float(mjd),max(mjd+10.)+float(mjd)]
nobs = n_elements(tobs)
tau=100.
if n_elements(tlag) eq 0 then  tlag= 300.
psi_width = 20.

dtpsi = psi_width/10.
tpsi_max = tlag + 5*psi_width
ngrid_psi = float(ceil(tpsi_max/dtpsi))
tpsi_true = dtpsi * findgen(ngrid_psi)
psi_true = exp(-(tpsi_true-tlag)^2/2./psi_width^2)
psi_true = psi_true/int_tabulated(tpsi_true,psi_true)


;Add a small amount of white noise.
noise_variance = .001^2
noise_vector = replicate(noise_variance,2*nobs)

;Now let us see if we can correctly back out the transfer function.
;Let's solve for psi in bins of width=w, ranging from t=0 to t=tmax
if ~keyword_set(w) then  w= 40.
tmax = 1000.
tmin = 0.
nbins = float(ceil((tmax-tmin)/w))
tpsi = [w/2+w*findgen(nbins)]


psi_in = exp(-(tpsi-tlag)^2/2./psi_width^2)
psi_in *= 0.
psi_in[0] = 1.0
;psi_in = interpol(psi,tobs[0:n_elements(psi)-1],tpsi)
psi_in = psi_in/int_tabulated(tpsi,psi_in)

;Next, compute the covariance matrices and the optimal estimator and all.



niter= 1000.
psi_avg = fltarr(nbins,niter)
for i=0L,niter-1 do begin
    print,'Iteration: ',string(i,form='(I0)')
    y = drw_sim_covar(tobs,tpsi_true,psi_true,sigma=sigma,mu=tau,$
                      noise=noise_vector,reset = (i eq 0))
    optimal_estimator,tobs,tpsi,y,w,dt,sigma=1.0,mu=tau,psi_in=psi_in,$
      psi_out=psi_out_intermediate,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,$
      noise=noise_vector
    optimal_estimator,tobs,tpsi,y,w,dt,sigma=1.0,mu=tau,$
      psi_in=psi_out_intermediate,psi_out=psi_out,$
      C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise_vector
    psi_avg[*,i] = -psi_out
;    plot,tpsi,total(psi_avg,2)/float(niter)
;    oplot,tpsi_true,psi_true,color=1.5e6

endfor
save,psi_avg,file=string('("estimator_test.mjd.tlag.",I05,".sav")',long(tlag))
prepare_plots,/color

filename = string(form='("../Plots/estimator_test.mjd.lownoise.tlag.",I05,".ps")',long(tlag))

psopen,filename,/color,xsize=6,ysize=6,/inches
plot,tpsi,total(psi_avg,2)/float(niter),xtitle='t (days)',ytitle='!7 w !6',$
  yr=[-1.0,1.5]
oplot,tobs,exp(-(tobs-tlag)^2/2./psi_width^2),color=200
;legend,['input','estimated'],line=0,color=[-1,200],/bottom,/right,box=0,charsize=2
psclose
prepare_plots,/reset



end

