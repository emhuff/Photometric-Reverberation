function C_CCcov,tobs,sigma,mu
;This function takes a list of observation times, and produces the
;covariance matrix of the continuum light curve.

nobs = n_elements(tobs)
ti = tobs # replicate(1,nobs)
tj = transpose(ti)

;By hypothesis, the covariance matrix of the continuum is given
;below. This is equation (2) from covariance.tex

C_cc = sigma^2*Exp(-(abs(ti-tj)/mu))

return,C_cc
end

function C_CLcov_m,tobs,tpsi,psi_m,sigma,mu,w,delta_t
;This function takes time indices and a guess for psi, along with the
;quasar timeseries parameters, and returns the m'th element of the
;continuum-line cross-covariance.
;
;To get the correct C_CL, sum this function over m. Note that psi and
;tpsi must both be scalars here.
if n_elements(tpsi) ne 1 then stop
if n_elements(psi_m) ne 1 then stop
nobs = n_elements(tobs)
ti = tobs # replicate(1,nobs)
tj = transpose(ti)
C_CL_m = fltarr(nobs,nobs)
dt  = ti-tj - tpsi

ind1 = where(dt gt w,ct1)
ind2 = where(dt ge -w AND dt le w,ct2)
ind3 = where(dt lt -w, ct3)

if ct1 gt 0 then begin
    y = dt[ind1]
    C_CL_m[ind1] = sigma^2*mu*psi_m/w* $
      exp(-w/mu - y/mu) * (exp(2*w/mu)-1)
endif

if ct2 gt 0 then begin
    y = dt[ind2]
    C_CL_m[ind2] = - sigma^2*mu*psi_m/w* $
      exp(-w/mu - y/mu) * $
      (1+ exp(2*y/mu) - 2*exp(w/mu+y/mu))
endif

if ct3 gt 0 then begin
    y = dt[ind3]
    C_CL_m[ind3] = sigma^2*mu*psi_m/w * $
      2*exp(y/mu)*sinh(w/mu)
endif

return,C_CL_m/2.
end


function C_LLcov_mn,tobs,tpsi_m,tpsi_n,psi_m,psi_n,sigma,mu,w,delta_t
if n_elements(tpsi_m) ne 1 then stop
if n_elements(psi_m) ne 1 then stop
if n_elements(tpsi_n) ne 1 then stop
if n_elements(psi_n) ne 1 then stop
nobs = n_elements(tobs)
ti = tobs # replicate(1,nobs)
tj = transpose(ti)

dt = ti-tj + (tpsi_m-tpsi_n)



C_LL_mn = fltarr(nobs,nobs)

ind1 = where(dt gt  w,ct1)
ind2 = where(dt lt -w,ct2)
ind3 = where(dt ge -w AND dt le 0, ct3)
ind4 = where(dt gt  0 AND dt le w, ct4)

if ct1 gt 0 then begin
    y = dt[ind1]
    C_LL_mn[ind1] = sigma^2*mu^2/w^2*psi_m*psi_n * $
      exp(-w/mu - y/mu)*(exp(w/mu)-1.)^2
endif

if ct2 gt 0 then begin
    y = dt[ind2]
    C_LL_mn[ind2] = sigma^2*mu^2/w^2*psi_m*psi_n * $
      exp(-w/mu + y/mu)*(exp(w/mu)-1.)^2
endif
if ct3 gt 0 then begin
    y = dt[ind3]
    C_LL_mn[ind3] = -sigma^2*mu/w^2*psi_m*psi_n * $
      exp(-w/mu-y/mu) * $
      (-mu - mu*exp(2*y/mu) + 2*mu*exp(w/mu+2*y/mu) - 2*w*exp(w/mu+y/mu) $
       - 2*y*exp(w/mu+y/mu))
endif
if ct4 gt 0 then begin
    y = dt[ind4]
    C_LL_mn[ind4] = sigma^2*mu/w^2*psi_m*psi_n * $
      exp(-w/mu-y/mu) * $
      (mu - 2*mu*exp(w/mu) + mu*exp(2*y/mu) + 2*w*exp(w/mu+y/mu) - $
       2*y*exp(w/mu+y/mu))
endif

return,C_LL_mn
end

pro optimal_estimator,tobs,tpsi,data,w,dt,sigma=sigma,mu=mu,$
                      psi_in=psi_in,psi_out=psi_out,$
                      C_CC=C_CC,C_CL=C_CL,C_LL=C_LL,$
                      noise=noise

nobs = n_elements(tobs)
npsi = n_elements(tpsi)

if ~keyword_set(psi_in) then psi_in = replicate(1./float(npsi),npsi)
if ~keyword_set(sigma) then sigma = 1.0
if ~keyword_set(mu) then mu = 100.

C_CL  = fltarr(nobs,nobs)
C_CLt = fltarr(nobs,nobs)
C_LL  = fltarr(nobs,nobs)


C_CC = C_CCcov(tobs,sigma,mu)
print,'making C_CL:'
for m = 0L,npsi-1 do begin
    print,string(form= '("Progress:",I02,"%")',float(m)/float(npsi)*100.)
    C_CL  += C_CLcov_m( tobs,tpsi[m],psi_in[m],sigma,mu,w,dt)
    C_CLt += C_CLcov_m(-tobs,tpsi[m],psi_in[m],sigma,mu,w,dt)
endfor
print,'Making C_LL:'
for m = 0L,npsi-1 do begin
    print,string(form= '("Progress:",I02,"%")',float(m)/float(npsi)*100.)
    for n=0L,npsi-1 do begin
        C_LL += C_LLcov_mn(tobs,tpsi[m],tpsi[n],psi_in[m],psi_in[n],sigma,mu,w,dt)
    endfor
endfor

Ncovar = diag_matrix(noise)

;Make the Big Covariance Matrix.

C = [[C_CC,C_CL],[transpose(C_CL),C_LL]]+Ncovar
Cinv = LA_invert(C,status=status,/double)
Fisher = dblarr(npsi,npsi)
q = dblarr(npsi)
f = dblarr(npsi)
;Make the Fisher matrix.
print,'Building the Fisher Matrix.'
for n = 0,npsi-1 do begin
    print,string(form= '("Progress:",I02,"%")',float(n)/float(npsi)*100.)
    for m = 0,npsi-1 do begin
        dCL_dpsi_m  = C_CLcov_m( tobs,tpsi[m],1.,sigma,mu,w,dt)
        dCL_dpsi_mt = C_CLcov_m(-tobs,tpsi[m],1.,sigma,mu,w,dt)
        dCL_dpsi_n  = C_CLcov_m( tobs,tpsi[n],1.,sigma,mu,w,dt)
        dCL_dpsi_nt = C_CLcov_m(-tobs,tpsi[n],1.,sigma,mu,w,dt)
        dLL_dpsi_m  = fltarr(nobs,nobs)
        dLL_dpsi_n  = fltarr(nobs,nobs)
        for i = 0,npsi-1 do begin
            dLL_dpsi_m += C_LLcov_mn(tobs,tpsi[m],tpsi[i],1.,psi_in[i],sigma,mu,w,dt)
            dLL_dpsi_n += C_LLcov_mn(tobs,tpsi[i],tpsi[n],psi_in[i],1.,sigma,mu,w,dt)
;function C_LLcov_mn,tobs,tpsi_m,tpsi_n,psi_m,psi_n,sigma,mu,w,delta_t
        endfor
;Now make dC_dm:
        dC_dm = [[C_CC*0.,dCL_dpsi_m],[dCL_dpsi_mt,dLL_dpsi_m]]
;And make dC_dn:
        dC_dn = [[C_CC*0.,dCL_dpsi_n],[dCL_dpsi_nt,dLL_dpsi_n]]
        Fisher[m,n] = 0.5*trace(Cinv # dC_dm # Cinv # dC_dn)
    endfor
    q[n] = 0.5*trace(Cinv # dC_dn # Cinv # Ncovar)
    f[n] = 0.5*transpose(data)# Cinv # dC_dn # Cinv # data
endfor

Finv = LA_invert(Fisher,status=status2,/double)
psi_out = Finv # (q-f)

end
