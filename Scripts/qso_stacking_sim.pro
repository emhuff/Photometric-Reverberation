pro qso_stacking_sim,nrealize=nrealize
if ~keyword_set(nrealize) then nrealize=1.


;Make a qso lightcurve. Start by making the grid nice and even.
ngrid = 1000.
tgrid = findgen(ngrid)
qso_lightcurve_sim,tgrid,tlag=100,c=c,l=l,transfer=transfer

;Next, we want to solve for psi. But we want to do it in
;relatively widely-spaced bins.
tpsi = findgen(10)*100.
data = [c,l]

optimal_estimator,tgrid,tpsi,data,psi=psi_out

stop

end
