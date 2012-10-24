function qso_varcat,nobj = nobj,rerun=rerun,input=input
if ~keyword_set(rerun) then rerun=301
qsocatfile = getenv('COADD_DATA_ROOT')+'/QSO_Reverb/Data/S82_quasars4EH.fits'
qsocat = mrdfits(qsocatfile,1)

if keyword_set(nobj) then begin
    if nobj ge n_elements(qsocat) then nobj = n_elements(qsocat)
    objind = random_indices(nobj,n_elements(qsocat))
    qsocat = qsocat[objind]
endif

runlist = sdss_runlist(rerun=301,/reload)
runlist = runlist[where(runlist.stripe eq 82)]
objlist = sdss_findallobj(ra,dec,run=runlist.run,$
                          parent=parent,child=child,rerun=301)

return,child
end
