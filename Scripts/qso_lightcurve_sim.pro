function fft_lightcurve,ngrid,sigma,tau
noise  = randomu(seed,ngrid)
fnoise = fft(noise)
k1 = findgen(ngrid/2+1)
k2 = -reverse(findgen(ngrid/2-1)+1)
kgrid = [temporary(k1),temporary(k2)]

pk = 4*sigma^2*tau/(1 + 2*!Pi*tau*kgrid^2)
pk = sigma^2/abs(kgrid)
pk[where(kgrid eq 0)] = sqrt(sigma)
phase = atan(fnoise,/phase)
pnoise = pk*exp(complex(0,1)*phase)
lightcurve = real_part(fft(pnoise,/inverse))

return,lightcurve
end

pro qso_lightcurve_sim,t_in,tau=tau,sigma=sigma,tlag=tlag,continuum=continuum,line=line,flux=flux,doplot=doplot,transfer_sigma = transfer_sigma,psi=transfer
if ~keyword_set(tau) then tau = 100.0
if ~keyword_set(sigma) then sigma = 1.0
if ~keyword_set(transfer_sigma) then transfer_sigma = 50.
if ~keyword_set(tlag) then tlag = 500.

t = t_in - min(t_in)

;First, find the minimum necessary grid spacing.
dt = t_in - shift(t_in,1)
dt = dt[1:(n_elements(dt)-1)]
nonil = where(dt ne 0,nilct)
if nilct eq 0 then stop
dt = dt[nonil]
dtmin = float(min(dt))

burn_interval = max([10*tau, 10*tlag+5*transfer_sigma])
interval = max(t_in)-min(t_in) + burn_interval
ngrid = float(ceil(interval/dtmin))
tgrid = dtmin*findgen(ngrid)
;The current value of the timeseries depends on previous values out to
;~tlag, so we build in a 'burn interval' that is long enough so that
;we don't accidentally change the noise properties.
lightcurve_oversampled = drw_sim(tgrid,tau,sigma)

;lightcurve_oversampled = fft_lightcurve(ngrid,sigma,tau)
continuum = interpol(lightcurve_oversampled,tgrid,t+burn_interval)

;Define a transfer function. Truncate so the convolution kernel isn't
;too big.
transfer = exp(-(tgrid-tlag)^2/2./(2.0*transfer_sigma))
transfer = transfer/total(transfer,/double)
transfer = transfer[where(tgrid lt tlag + 5*transfer_Sigma)]

line_oversampled = convol(lightcurve_oversampled,transfer,/edge_truncate,center=0)
line = interpol(line_oversampled,tgrid,t+burn_interval)

if keyword_set(Doplot) then begin
;    psopen,'sample_lightcurve',xsize=6,ysize=6,/inches,/color,/enc
;    prepare_plots,/color
;    device,/color
    loadct,6
    plot,tgrid,lightcurve_oversampled,xtitle='time (arbitrary)'
    oplot,tgrid,line_oversampled,color=150
;    psclose
endif


end
