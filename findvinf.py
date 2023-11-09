import findtca
import numpy as np
import spiceypy as sp

class FINDVINF:
  def __init__(self
              ,tcaInp
              ,VinfTargetArg=None
              ,VinfEtOffsArg=None
              ,VinfEtStepArg=None
              ,VinfEtNumArg=None
              ,etRtnOffset=0e0
              ,kernels=[]
            ):
    self.success, self.status, self.msg = False, 'Vinf CALC SKIPPED', ''
    if VinfTargetArg is None: return

    etOffDefault = sp.spd() * 3.25

    self.VinfTarget = str(VinfTargetArg).strip()

    try   : self.VinfEtOffs = float(VinfEtOffsArg)
    except: self.VinfEtOffs = etoffDefault

    try   : self.VinfEtStep = float(VinfEtStepArg)
    except: self.VinfEtStep = 600e0

    try   : self.VinfEtNum = int(VinfEtNumArg)
    except: self.VinfEtNum = 24

    self.status = 'Vinf CALC FAILED'

    try:
      findtca.ftca_furnsh(kernels)

      caOut=findtca.findtca(self.VinfTarget,tcaInp.Obs
                           ,utcEst='JD'+str(tcaInp.JDTdbTca,2).strip()
                           )

      etTca = caOut.et
      spd = sp.spd()
      n1 = self.VinfEtNum
      iw1 = np.arange(n1,dtype=np.int32)
      n = n1*3L
      iw23 = n+np.arange(n1*2L,dtype=np.int32)
      t0s =  -rtn.VinfEtOffs + np.arange(n)*36d2
      ts = etTca+t0s
      vsTarg = np.zeros(3*n,dtype=np.float64).reshape((-1,3,))
      vsObs = np.zeros(3*n,dtype=np.float64).reshape((-1,3,))
      try   : targNum = sp.bodn2c(caOut.Target)
      except: targNum = int(caOut.Target)
      if (targNum >= 100 and targNum < 1000
         ) or (targNum >= 248695800 and targNum < 248695900
         ):
        baryID = targNum // 100
      else:
        baryID = targNum

      baryName = sp.bodc2s(baryID)

      for i in range(n):
        vsTarg[i] = sp.spkezr(caOut.Target,ts[i],'j2000','none',baryName)
        vsObs[i] = sp.spkezr(caOut.Obs,ts[i],'j2000','none',baryName)
      vs = vsTarg - vsObs
      ABs = dindgen(2,3)
ABsObs = dindgen(2,3)
for i = 0L,2L do begin
  vis = vs[i,iw1]
  ABs[*,i] = linfit(t0s[iw1],vis[*],/double,yfit=yfits)
  vis = vsObs[i,iw1]
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
;;; +Zbpl in J2000 frame  = > TOF   => Obs asymptotic velocity wrt Targ
;;; +Ybpl in J2000 frame  = > BdotT => [Target position wrt Obs] cross +Zbpl
;;; +Xbpl in J2000 frame  = > BdotR => Obs position vec wrt Targ @ TCA

cspice_vhat,  velObsTca, uvTof                 ;;; TOF    = > +Zbpl
cspice_ucrss,  pos, uvTof, uvBdotT             ;;; BdotT  = > +Ybpl
cspice_ucrss,  uvBdotT, uvTof, uvBdotR         ;;; BdotR  = > +Xbpl

;;; [+Xbpl, +Ybpl, +Zbpl]  = > Matrix to convert J2000 vectors to B-plane frame

mtx_j2b = double( [uvBdotR, uvBdotT, uvTof], 0, 3, 3)

cspice_recrad,uvTof,rTof,raTof,decTof
cspice_recrad,uvBdotR,rBdotR,raBdotR,decBdotR

dpr = cspice_dpr()

raTof = raTof * dpr
decTof = decTof * dpr

raBdotR = raBdotR * dpr
decBdotR = decBdotR * dpr

;;; Xform from barycentric target RTN frame to J2000
;;; R  = > Radial from Barycenter to target
;;; N  = > Normal to obit plane, R cross velocity
;;; T  = > Along track, N cross R

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
message,'OK',continue = debug

    except:
       ftca_furnsh(kernels, /unload)
       self.msg = traceback.format_exc()

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FVINF_TEST:  FINDTCA support function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;n
function fvinf_test, targ, obs, kernels = kernels $
             , debug = debug

  ftca_furnsh, kernels, debug = debug

  TcaOut = findtca(targ, obs, debug=debug)

  ftca_furnsh, kernels, /unload, debug = debug

  return, TcaOut

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ks = ['nh_ref_v_od059b_encounter_only.bsp']

!quiet = 1
testEr = 0L
catch,testEr

if testEr eq 0L then begin
  pcaOut = fvinf_test('999',kernels=ks)
  pvinfOut = findvinf(pcaOut,'999',kernels=ks,/debug,etRtnOffset=-6d2)
  cvinfOut = findvinf(pcaOut,'901',kernels=ks,/debug)
  nvinfOut = findvinf(pcaOut,'902',kernels=ks,/debug)
  hvinfOut = findvinf(pcaOut,'903',kernels=ks,/debug)
  catch,/cancel
endif else begin
  catch,/cancel
  message,/continue,!error_state.msg
  message,/reset
endelse
help,/st,pvinfOut

end
