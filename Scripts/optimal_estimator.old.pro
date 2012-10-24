pro optimal_estimator,tobs,tpsi,data,sigma=sigma,w=w,s0=s0,psi=psi
;sigma is the typical error on the observed fluxes.
if ~keyword_Set(sigma) then sigma=1.
if ~keyword_set(w)  then w=1.
if ~keyword_set(s0) then s0=1.
nobs = n_elements(tobs)
npsi = n_elements(tpsi)
;Assume a power spectrum.
;P(w) = S0/( 1 + (S1/w)^2 )
;For this, then, xi(ti,tj) = xi(ti-tj) = exp(-abs(ti-tj)/w)
ti = tobs # replicate(1,nobs)
tj = transpose(ti)
tij = ti-tj                               

tm  = tpsi # replicate(1,npsi)
tn  = transpose(tm)
tmn = tm-tn

psi_guess = replicate(1.,npsi)

xi = s0*exp(-(abs(tij)/w))

;Compute the covariance matrix between the observation times.
cov_CC = sigma^2*xi;
dcov_CC_dpsi = fltarr(nobs,nobs)
cov_CL = fltarr(nobs,nobs)
dcov_CL_dpsi = fltarr(nobs,nobs)
cov_LL = fltarr(nobs,nobs)
dcov_LL_dpsi = fltarr(nobs,nobs)
for i=0,nobs-1 do begin
    print,i
    for j=0,nobs-1 do begin
        cov_CL[i,j] = total(psi_guess*s0*exp(-(abs(tpsi-tij[i,j]))/w))
        dcov_CL_dpsi[i,j] = total(s0*exp(-(abs(tpsi-tij[i,j]))/w))
        for m=0,npsi-1 do begin
            for n=0,npsi-1 do begin
                cov_LL[i,j] += psi_guess[m]*psi_guess[n]*s0*exp(-(abs(tmn[m,n]-tij[i,j])/w))
                dcov_LL_dpsi[i,j] += psi_guess[n]*s0*exp(-(abs(tmn[m,n]-tij[i,j])/w))
            endfor
        endfor
    endfor
endfor


;Assemble the disparate covariance matrices into the full data
;covariance matrix.
cov   = [[[cov_CC],[cov_CL]],transpose([cov_CL,cov_LL])]
noise = diag_matrix(replicate(sigma,2*nobs))
cov = cov+noise
dcovm  = [[[dcov_CC_dpsi],[dcov_CL_dpsi]],transpose([dcov_CL_dpsi,dcov_LL_dpsi])]
dcovn  = [[[dcov_CC_dpsi],[dcov_CL_dpsi]],transpose([transpose(dcov_CL_dpsi),transpose(dcov_LL_dpsi)])] 
cinv  = LA_invert(cov,status=status)

if status ne 0 then stop

;Compute the pieces of the estimator.

Fisher = fltarr(npsi,npsi);
q = fltarr(npsi)
f = fltarr(npsi)
for m = 0,npsi-1 do begin
    for n=0,npsi-1 do begin
        Fisher[m,n] = trace(cinv # dcovm # cinv # dcovn)/2.
    endfor
endfor
f = trace(cinv # dcovm # cinv # noise)/2.
q = transpose(data) # (cinv # dcovm # cinv) # data/2.

stop
;Finally, assemble the Estimate!
psi = fltarr(npsi)
Finv = LA_invert(Fisher,status=status2)
if status2 gt 0 then stop

for m=0,npsi-1 do begin
    psi[m] = total(Finv # (q-f),2)
endfor

end
