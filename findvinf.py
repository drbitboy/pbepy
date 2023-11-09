def findVInf(tcaInp
            ,VinfTargetArg $
            ,VinfEtOffsArg $
            ,VinfEtStepArg $
            ,VinfEtNumArg $
            ,etRtnOffset=etRtnOffsetArg $
            ,kernels=kernels $
            ):

dodebug = keyword_set(debug)
nodebug = dodebug ? 0b : 1b

resolve_all,/continue

rtn = { success: 0b, status: 'Vinf CALC SKIPPED', msg:'' }

if n_elements(VinfTargetArg) ne 1L then return, rtn

etRtnOffset = n_elements(etRtnOffsetArg) eq 1L ? double(etRtnOffsetArg[0]) : 0d0

etOffDefault = cspice_spd() * 3.25

rtn = create_struct( rtn $
, 'VinfTarget', strtrim(VinfTargetArg[0],2) $
, 'VinfEtOffs' $
, n_elements(VinfEtOffsArg) eq 1 ? double(VinfEtOffsArg[0]) : etOffDefault $
, 'VinfEtStep' $
,  n_elements(VinfEtStepArg) eq 1 ? double(VinfEtStepArg[0]) : 600d0 $
, 'VinfEtNum' $
, n_elements(VinfEtNumArg) eq 1 ? long(VinfEtNumArg[0]) : 24 $
)

rtn.status = 'Vinf CALC FAILED'

testEr=0L
catch,testEr

if testEr ne 0L then begin
  catch,/cancel
  ftca_furnsh, kernels, /unload, debug=debug
  rtn.msg=!error_state.msg
  message,/reset
  return, rtn
endif

ftca_furnsh, kernels, debug=debug

caOut=findtca(rtn.VinfTarget,tcaInp.Obs $
             ,utcEst='JD'+strtrim(tcaInp.JDTdbTca,2))

etTca=caOut.et
spd=cspice_spd()
n1=rtn.VinfEtNum
iw1=lindgen(n1)
n=n1*3L
iw23=n+lindgen(n1*2L)
t0s= -rtn.VinfEtOffs + dindgen(n)*36d2
ts=etTca+t0s
vsTarg=dblarr(3,n)
vsObs=dblarr(3,n)
cspice_bodn2c,caOut.Target,targNum,f
if not f then targNum = long(caOut.Target)
baryID = (targNum ge 100 and targNum lt 1000) $
      OR (targNum ge 248695800 and targNum lt 248695900) $
       ? (targNum / 100L) $
       : targNum
cspice_bodc2s,baryID,baryName

for i=0L,n-1L do begin
  cspice_spkezr,caOut.Target,ts[i],'j2000','none',baryName,st,ltim
  vsTarg[*,i]=st[0:2]
  cspice_spkezr,caOut.Obs,ts[i],'j2000','none',baryName,st,ltim
  vsObs[*,i]=st[0:2]
endfor
vs = vsTarg - vsObs
ABs = dindgen(2,3)
ABsObs = dindgen(2,3)
for i=0L,2L do begin
  vis=vs[i,iw1]
  ABs[*,i] = linfit(t0s[iw1],vis[*],/double,yfit=yfits)
  vis=vsObs[i,iw1]
  ABsObs[*,i] = linfit(t0s[iw1],vis[*],/double,yfit=yfits)
endfor

vel = -transpose(ABs[1,*])             ;;; Target-relative Obs velocity
velObs = transpose(ABsObs[1,*])        ;;; Bary-relative Obs velocity
avgEt = mean(ts,/double)               ;;; Mean Time
avgpos = total(vs,2,/double)/n         ;;; Mean S/C-relative target position
avgPosObs = total(vsObs,2,/double)/n   ;;; Mean Bary-relative Obs position

;;; find B-plane TCA 'BTCA' i.e. when B-plane passes through target

et = etTca
dt = 0d0
ltim2 = 0d0
iter = 0L
while ( iter lt 20L and abs(dt) gt 1d-3 ) or iter eq 0L do begin
  iter = iter + 1
  et = et + dt
  cspice_spkezr,caOut.Target,et,'j2000','none',baryName,st,ltim ;;;Target @ BTCA
  targWrtBary = st[0:2]
  velTargWrtBary = st[3:5]
  obsWrtBary = avgPosObs + velObs * (et + ltim2 - avgEt)
  pos = targWrtBary - obsWrtBary

  ;;; do one iteration on light time
  ltim2 = cspice_vnorm(pos) / cspice_clight()
  obsWrtBary = avgPosObs + velObs * (et + ltim2 - avgEt)
  pos = targWrtBary - obsWrtBary

  vel = velObs - velTargWrtBary
  ;;; N.B. This is Vinfinity, not the observer@TCA velocity wrt the Target@TCA
  dt = cspice_vdot(vel,pos) / cspice_vdot(vel,vel)  ;;; delta ET to BTCA
endwhile
if dodebug then help,iter,dt,ltim2,et-etTca

velObsTca = vel

;;; Xform from J2000 to B-plane frame, subscript bpl
;;;
;;; +Zbpl in J2000 frame => TOF   => Obs asymptotic velocity wrt Targ
;;; +Ybpl in J2000 frame => BdotT => [Target position wrt Obs] cross +Zbpl
;;; +Xbpl in J2000 frame => BdotR => Obs position vec wrt Targ @ TCA

cspice_vhat,  velObsTca, uvTof                 ;;; TOF   => +Zbpl
cspice_ucrss,  pos, uvTof, uvBdotT             ;;; BdotT => +Ybpl
cspice_ucrss,  uvBdotT, uvTof, uvBdotR         ;;; BdotR => +Xbpl

;;; [+Xbpl, +Ybpl, +Zbpl] => Matrix to convert J2000 vectors to B-plane frame

mtx_j2b = double( [uvBdotR, uvBdotT, uvTof], 0, 3, 3)

cspice_recrad,uvTof,rTof,raTof,decTof
cspice_recrad,uvBdotR,rBdotR,raBdotR,decBdotR

dpr = cspice_dpr()

raTof = raTof * dpr
decTof = decTof * dpr

raBdotR = raBdotR * dpr
decBdotR = decBdotR * dpr

;;; Xform from barycentric target RTN frame to J2000
;;; R => Radial from Barycenter to target
;;; N => Normal to obit plane, R cross velocity
;;; T => Along track, N cross R

etRtn = tcaInp.et + etRtnOffset
cspice_etcal,etRtn,rtnTDBCaldate
cspice_spkezr,caOut.Target,etRtn,'j2000','none',baryName,stBary2TargTca,ltim

posBary2Targ = stBary2TargTca[0:2]
velBary2Targ = stBary2TargTca[3:5]

cspice_vhat, posBary2Targ, uvRb2tt
cspice_ucrss, uvRb2tt, velBary2Targ, uvNb2tt
cspice_ucrss, uvNb2tt, uvRb2tt, uvTb2tt

rtn = create_struct( rtn $
, 'mtx_j2b', mtx_j2b $
, 'raTof', raTof $
, 'decTof', decTof $
, 'raBdotR', raBdotR $
, 'decBdotR', decBdotR $
, 'bodyTcaBary2Targ', tcaInp.Target $
, 'tdbBary2Targ', tcaInp.CALTIMESPICE $
, 'posBary2Targ', posBary2Targ $
, 'velBary2Targ', velBary2Targ $
, 'mtx_RTNtarg2j', transpose( double([ uvRb2tt, uvTb2tt, uvNb2tt], 0, 3, 3)) $
, 'tdbRTNCaldate', rtnTDBCaldate + ' TDB (SPICE ETCAL)' $
, 'tdbRTN_sPastJ2k', etRtn $
, 'Rtn_Offset_s', etRtnOffset $
)

rtn.success = 1b
rtn.status = 'Vinf CALC SUCCEEDED'
message,'OK',continue=debug

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FVINF_TEST:  FINDTCA support function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;n
function fvinf_test, targ, obs, kernels=kernels $
             , debug=debug

  ftca_furnsh, kernels, debug=debug

  TcaOut = findtca(targ, obs, debug=debug)

  ftca_furnsh, kernels, /unload, debug=debug

  return, TcaOut

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ks=['nh_ref_v_od059b_encounter_only.bsp']

!quiet=1
testEr=0L
catch,testEr

if testEr eq 0L then begin
  pcaOut=fvinf_test('999',kernels=ks)
  pvinfOut=findvinf(pcaOut,'999',kernels=ks,/debug,etRtnOffset=-6d2)
  cvinfOut=findvinf(pcaOut,'901',kernels=ks,/debug)
  nvinfOut=findvinf(pcaOut,'902',kernels=ks,/debug)
  hvinfOut=findvinf(pcaOut,'903',kernels=ks,/debug)
  catch,/cancel
endif else begin
  catch,/cancel
  message,/continue,!error_state.msg
  message,/reset
endelse
help,/st,pvinfOut

end
