function drw_sim_covar,tobs,tpsi,psi,$
                       sigma=sigma,mu=mu,reset=reset,noise=noise
COMMON MASTER_COVARIANCE,M_master
if  keyword_set(reset) then undefine,m_master
if ~keyword_set(sigma) then sigma=1.
if ~keyword_set(mu) then mu = 100.

if keyword_set(M_master) then M = M_master $
else begin
    nobs = n_elements(tobs)
    w = tpsi[1]-tpsi[0]
    if ~keyword_set(noise) then noise = replicate(0.01,2*nobs)
    
    reverb_covariance,tobs,tpsi,psi,w,$
      sigma=sigma,covar=Cov,$
      noise=noise,C_CL=C_CL
    

    ti = [findgen(2*nobs)] # replicate(1.,2*nobs)
    tj = transpose(ti)
    M = cov
    mask = float(ti ge tj)
    la_choldc,M,status=status,/upper,/double
;    if status ne 0 then stop
    M = M * mask
    M_master = M
endelse

dev = randomn(seed,2*n_elements(tobs))
y = M # dev


return,y
end

