function drw_sim,tobs,tau,sigma, M=M
nobs = n_elements(tobs)
ti = tobs # replicate(1.,nobs)
tj = transpose(ti)
cov = sigma^2 *exp(-abs(ti-tj)/tau)
M = cov
mask = ti ge tj
la_choldc,M,status=status,/upper,/double
M = M * mask
dev = randomn(seed,nobs)
if status ne 0 then stop
y = M # dev
return,y
end
