pro s82_simulate,tlag=tlag,w=w

;Generate light curves using the actual Stripe 82 cadences.
restore,'sdss_S82_cadences.sav'
nobs_vec = total(mjd_obs ge 0,2)
usable = where(nobs_vec gt 20 AND zlist gt 0)
zlist = zlist[usable]
mjd_obs = mjd_obs[usable,*]

;Note that this contains two arrays:
; 1. mjd_obs[i,j] is the mjd of the j'th observation of the i'th
;    quasar. It is equal to -1 after the last observation.
; 2. zlist is an array with redshifts randomly sampled from the list
;    of spectroscopic observations on S82
;
; This code generates a new realization of a quasar lightcurve at the
; given redshift. Each realization has a different observation cadence.


tau = 100.
if n_elements(tlag) eq 0 then  tlag= 150.
psi_width = 60.

dtpsi = psi_width/10.
tpsi_max = tlag + 5*psi_width
ngrid_psi = float(ceil(tpsi_max/dtpsi))
tpsi_true = dtpsi * findgen(ngrid_psi)
psi_true = exp(-(tpsi_true-tlag)^2/2./psi_width^2)
psi_true = 100.0*psi_true/int_tabulated(tpsi_true,psi_true)

;Now let us see if we can correctly back out the transfer function.
;Let's solve for psi in bins of width=w, ranging from t=0 to t=tmax
if ~keyword_set(w) then  w= 30.
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



niter= 100.
qsoind = random_indices(n_elements(zlist),niter)

psi_avg = fltarr(nbins,niter)
condrecord = lonarr(3,niter)
for i=0L,niter-1 do begin

;   Choose the observation time and the redshift.
    tobs = mjd_obs[qsoind[i],*]
    tobs = tobs[where(tobs ge 0)]
    tobs = tobs - min(tobs)
    z = zlist[qsoind[i]]

    print,'Iteration: ',string(i,form='(I0)'),'  nobs: ',string(n_elements(tobs),form='(I0)')
;Add a realistic amount of white noise 
;(the median i-band SNR for these quasars is ~20)
    noise_variance = .01^2
    noise_vector = replicate(noise_variance,2*n_elements(tobs))
    
    y = drw_sim_covar(tobs/(1+z),tpsi_true,psi_true,sigma=sigma,mu=tau,$
                      noise=noise_vector,/reset,cond=cond1)
    optimal_estimator,tobs/(1+z),tpsi,y,w,dt,sigma=1.0,mu=tau,psi_in=psi_in,$
      psi_out=psi_out_intermediate,C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,$
      noise=noise_vector,cond=cond2
    optimal_estimator,tobs/(1+z),tpsi,y,w,dt,sigma=1.0,mu=tau,$
      psi_in=psi_out_intermediate,psi_out=psi_out,$
      C_CC=CCcov,C_CL=CLcov,C_LL=LLcov,noise=noise_vector,cond=cond3
    psi_avg[*,i] = -psi_out
    plot,tpsi,-psi_out;total(psi_avg,2)/float(niter)
    oplot,tpsi_true,psi_true,color=1.5e6
    xyouts,0.8,0.8,cond1+cond2+cond3,/norm,charsize=1.5
    print,'Matrix conditioning: ',cond1,cond2,cond3
    condrecord[*,i] = [cond1,cond2,cond3]
endfor
save,tpsi,psi_avg,file=string(form='("sdss_sim.tlag.",I05,".niter.",I0,".sav")',long(tlag),long(niter))
prepare_plots,/color
filename = string(form='("../Plots/sdss_sim.tlag.",I05,".niter.",I0,".ps")',long(tlag),long(niter))
stop
psopen,filename,/color,xsize=6,ysize=6,/inches
plot,tpsi,total(psi_avg,2)/float(niter),xtitle='t (days)',ytitle='!7 w !6',$
  yr=[-1.0,1.5]
oplot,tobs,exp(-(tobs-tlag)^2/2./psi_width^2),color=200
;legend,['input','estimated'],line=0,color=[-1,200],/bottom,/right,box=0,charsize=2
psclose
prepare_plots,/reset



end

