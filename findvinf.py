import findtca
import numpy as np
import spiceypy as sp

try: clight
except:
  clight = sp.clight()
  spd = sp.spd()
  dpr = sp.dpr()

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

    etOffDefault = spd * 3.25

    self.VinfTarget = str(VinfTargetArg).strip()

    try   : self.VinfEtOffs = float(VinfEtOffsArg)
    except: self.VinfEtOffs = etOffDefault

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
      n1 = self.VinfEtNum
      iw1 = np.arange(n1,dtype=np.int32)
      n = n1*3
      iw23 = n+np.arange(n1*2,dtype=np.int32)
      ts = (np.arange(n)*36e2) + (etTca - rtn.VinfEtOffs)
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

      vsTarg = np.vstack(
         [sp.spkezr(caOut.Target,t,'j2000','none',baryName)[0][:2]
          for t in ts
         ])
      vsObs = np.vstack(
         [sp.spkezr(caOut.Obs,t,'j2000','none',baryName)[0][:2]
          for t in ts
         ])

      vs = vsTarg - vsObs

      a = np.vstack((np.ones(n),ts,))

      ### Linear fit time vs. position => Target-relative Obs velocity
      ### Linear fit time vs. position => Bary-relative Obs vel => vInf
      ### Mean Time
      ### Mean Bary-relative Obs position
      vel = - np.linalg.lstsq(a,vs,rcond=None)[0][1]
      vInf = np.linalg.lstsq(a,vsObs,rcond=None)[0][1]
      avgEt = ts.mean()
      avgPosObs = vsObs.mean(axis=0)

      ### find B-plane TCA 'BTCA' i.e. when target is in B-plane

      et = etTca
      dt = 0e0
      ltim2 = 0e0
      iter = 0
      while ( iter < 20 and abs(dt) > 1e-3 ) or iter == 0:
        iter = iter + 1
        et = et + dt
        ###Target @ BTCA
        ltim,st = sp.spkezr(caOut.Target,et,'j2000','none',baryName)
        targWrtBary = st[0:2]
        velTargWrtBary = st[3:5]
        obsWrtBary = avgPosObs + vInf * (et + ltim2 - avgEt)
        pos = targWrtBary - obsWrtBary

        ### do one iteration on light time
        ltim2 = sp.vnorm(pos) / clight
        obsWrtBary = avgPosObs + vInf * (et + ltim2 - avgEt)
        vTarg = targWrtBary - obsWrtBary

        vel = vInf - velTargWrtBary
        ### N.B. This is Vinfinity, not the observer@TCA velocity wrt the Target@TCA
        dt = sp.vdot(vel,vTarg) / sp.vdot(vel,vel)  ### delta ET to BTCA

      vInfTca = vel

      ### Xform from J2000 to B-plane frame, subscript bpl
      ###
      ### +Zbpl in J2000 frame  = > TOF   => Obs asymptotic velocity wrt Targ
      ### +Ybpl in J2000 frame  = > BdotT => [Target position wrt Obs] cross +Zbpl
      ### +Xbpl in J2000 frame  = > BdotR => Obs position vec wrt Targ @ TCA

      uvTof = sp.vhat(vInfTca)                      ### TOF    = > +Zbpl
      uvBdotT = sp.ucrss(vTarg, uvTof)              ### BdotT  = > +Ybpl
      uvBdotR = sp.ucrss(uvBdotT, uvTof)            ### BdotR  = > +Xbpl

      ### [+Xbpl, +Ybpl, +Zbpl]  = > Matrix to convert J2000 vectors to B-plane frame

      mtx_j2b = numpy.vstack((uvBdotR, uvBdotT, uvTof,))

      rTof,raTof,decTof = sp.recrad(uvTof)
      rBdotR,raBdotR,decBdotR - sp.recrad(uvBdotR)


      raTof = raTof * dpr
      decTof = decTof * dpr

      raBdotR = raBdotR * dpr
      decBdotR = decBdotR * dpr

      ### Xform from barycentric target RTN frame to J2000
      ### R  = > Radial from Barycenter to target
      ### N  = > Normal to orbit plane, R cross velocity
      ### T  = > Along track, N cross R

      etRtn = tcaInp.et + etRtnOffset
      rtnTDBCaldate = sp.etcal(etRtn)
      stBary2TargTca = sp.spkezr(caOut.Target,etRtn,'j2000','none',baryName)[0]

      posBary2Targ = stBary2TargTca[:3]
      velBary2Targ = stBary2TargTca[3:]

      uvRb2tt = sp.vhat(posBary2Targ)
      uvNb2tt = sp.ucrss(uvRb2tt, velBary2Targ)
      uvTb2tt = sp.ucrss(vNb2tt, uvRb2tt)

      self.mtx_j2b = mtx_j2b
      self.raTof = raTof
      self.decTof = decTof
      self.raBdotR = raBdotR
      self.decBdotR = decBdotR
      self.bodyTcaBary2Targ = tcaInp.Target
      self.tdbBary2Targ = tcaInp.CALTIMESPICE
      self.posBary2Targ = posBary2Targ
      self.velBary2Targ = velBary2Targ
      self.mtx_RTNtarg2j = numpy.vstack((uvRb2tt, uvTb2tt, uvNb2tt,)).T
      self.tdbRTNCaldate = rtnTDBCaldate + ' TDB (SPICE ETCAL)'
      self.tdbRTN_sPastJ2k = etRtn
      self.Rtn_Offset_s = etRtnOffset

      self.success = True
      self.status = 'Vinf CALC SUCCEEDED'
      self.msg = 'OK'

    except:
       ftca_furnsh(kernels, unload=True)
       self.msg = traceback.format_exc()

"""
########################################################################
### FVINF_TEST:  FINDTCA support function
#######################################################################n
function fvinf_test, targ, obs, kernels = kernels $
             , debug = debug

  ftca_furnsh, kernels, debug = debug

  TcaOut = findtca(targ, obs, debug=debug)

  ftca_furnsh, kernels, /unload, debug = debug

  return, TcaOut

end

########################################################################

ks = ['nh_ref_v_od059b_encounter_only.bsp']

!quiet = 1
testEr = 0
catch,testEr

if testEr eq 0 then begin
  pcaOut = fvinf_test('999',kernels=ks)
  pvinfOut = findvinf(pcaOut,'999',kernels=ks,/debug,etRtnOffset=-6e2)
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
"""
