;function drw_sim,ngrid,dt,tau,sigma,b
;;This is supposed to generate a damped random walk timeseries.
;;
;;tau   -- this is the timescale over which the influence of fluctuations
;;         remain
;;sigma -- this is the variance of the series
;;b     -- the mean of the series.
;
;
;;We can test this by calculating the autocorrelation, and comparing to
;;e^(-t/tau). It works if ngrid*dt >> tau.
;
;drw = fltarr(ngrid)
;
;db = dt*randomn(seed,ngrid)
;tgrid = findgen(ngrid)*dt
;drw = b*tau*(1-exp(-tgrid/tau))
;for i =1L,ngrid-1 do begin
;    t = i*dt
;    drw[i]  +=sigma * total(exp(-(t-tgrid[0:i])/tau)*db[0:i])
;endfor
;
;return,drw
;end


function drw_sim,tobs,tau,sigma
nobs = n_elements(tobs)
ti = tobs # replicate(1.,nobs)
tj = transpose(ti)
cov = sigma^2 *exp(-abs(ti-tj)/tau)
M = cov
LA_choldc,M,status=status,/upper
dev = randomn(seed,nobs)
if status ne 0 then stop
y = M # dev
return,y
end
