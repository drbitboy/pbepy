import sys
import numpy
import datetime
import traceback
import collections
import spiceypy as sp

import findvinf

werr = sys.stderr.write
try: clight
except:
  clight = sp.clight()
  spd = sp.spd()


#######################################################################n
def ftca_furnsh(kernels, unload=False):
  """
ftca_furnsh:  FINDTCA support procedure,
              load/unload multiple SPICE kernels

"""
  func = unload and sp.unload or sp.furnsh
  for kernel in kernels:
    try:
      func(kernel)
    except:
      traceback.print_exc()


#######################################################################n
def ftca_getrange(targ, obs, et, ltcorrIn='LT',debug=False):
  """
ftca_getrange:  FINDTCA support function, get range from observer (obs)
                to target (targ) at et (TDB) using light time correction
                (ltcorrIn=).  Also returns state vector (float[6]) from
                the observer to target, light time, light time corr used

"""
  try:
    stOut, ltime = sp.spkezr(targ,et,'J2000',ltcorrIn,obs)
    return ltime*sp.clight(),stOut,ltime,ltcorrIn
  except:
    if debug: traceback.print_exc()

  return -1e0,numpy.array([0e0]*6),-1e0/clight,ltcorrIn


########################################################################
def ftca_getopt(defawlt, arg, poolvars):
  """
FTCA_GETOPT:  FINDTCA support function, select option
              via argument choices in decreasing priority:
              ARG       Function argument
              POOLVARS  SPICE kernel pool variable(s)
              DEFAWLT   Default value

              return value,which

"""

  if not (arg is None): return arg,'ARGUMENT'

  if not isinstance(poolvars,collections.Sequence):
    return ftcs_getopt(defawlt,arg,[poolvars])

  try: typ = defawlt.dtype.name
  except:
    typ=((not (type(defawlt) is str))
         and isinstance(defawlt,collections.Sequence)
         and type(defawlt[0])
         or type(defawlt)
        ).__name__

  for poolvar in poolvars:
    try:
      if typ.startswith('str')    : Vals = sp.gcpool(poolvar,0,999,999)
      elif typ.startswith('float'): Vals = sp.gdpool(poolvar,0,999)
      else                        : Vals = sp.gipool(poolvar,0,999)

      if len(Vals) == 1: Vals = Vals.pop()
      return Vals,'KERNELPOOL'
    except:
      continue

  return defawlt,'DEFAULT'


########################################################################
def julday():
  """Approximate Julian of right now"""
  return ((datetime.datetime.now()-datetime.datetime(2000,1,1,12)
          ).total_seconds()/spd)+sp.j2000()


########################################################################
########################################################################
### FINDTCA:  Main function
########################################################################
########################################################################
class FINDTCA:
  def __init__(self
              ,targArg=None
              ,obsArg=None
              ,utcEstArg=None
              ,VinfTarg=None
              ,VinfEtOffset=None
              ,VinfEtStep=None
              ,VinfEtNum=None
              ,ltcorrIn='LT'
              ,kernels=[]
              ,nounloadkernels=False
              ):

    ftca_furnsh(kernels)

    self.Target = targ = ftca_getopt('PLUTO', targArg, 'FINDTCA_TARGET')[0]
    self.Obs = obs = ftca_getopt('-98', obsArg, 'FINDTCA_OBSERVER')[0]
    self.lt = ltcorrIn

    ### use UTC estimate if given as argument,
    ### else use kernel pool variable FINDTCA_UTCEST
    ### else use NOW

    defawlt = f'JD {julday():.1f} TDB'
    utcEst,wUtcEst = ftca_getopt(defawlt, utcEstArg, 'FINDTCA_UTCEST')

    ### Convert delEst to s

    delEst = spd * (wUtcEst == 'DEFAULT' and 7e0 or 1e0)

    try:

      ### Find out if LEAPSECONDS have been loaded
      ### Convert UTC estimate to ET (s past J2000)

      try:
        etEst = sp.str2et(utcEst)
      except:
        etEst,tpErrmsg = sp.tparse(utcEst)
        werr('WARNING:  NO LEAPSECOND KERNEL AVAILABLE\n')
        if tpErrmsg.strip() != '': werr(f'TPARSE:  {tpErrmsg}\n')

      ### Go up to 1000 steps either side of the UTC estimate, find the
      ### first three in a row when minimum range is positive and span
      ### closest approach

      iCenter = 1000
      ranges = -numpy.ones(iCenter * 2 + 1)

      iMin = -1
      i = 0
      while i <= iCenter and iMin < 0:

        dt = delEst * i

        ### Extend calculated ranges forward in time
        i0 = iCenter+i
        ranges[i0] = ftca_getrange(targ,obs,etEst+dt,ltcorrIn=self.lt)[0]
        #print({i0:(ranges[i0],targ,obs,sp.etcal(etEst+dt),)})
        if ranges[i0] > 4e7: ranges[i0] = -1e0

        ### If three ranges backward from i0 have positive ranges, ...
        if min(ranges[i0-2:i0+1]) > 0e0:
          ### ... and middle range is <= end ranges ...
          if min(ranges[i0-2:i0+1:2]) >= ranges[i0-1]:
            iMin = i0 - 1                       ### ... use middle range
            continue                            ### skip backward step

        ### Extend calculated ranges backward in time
        i0 = iCenter-i
        ranges[i0] = ftca_getrange(targ,obs,etEst-dt,ltcorrIn=self.lt)[0]
        #print({i0:(ranges[i0],targ,obs,sp.etcal(etEst+dt),)})
        if ranges[i0] > 4e7: ranges[i0] = -1e0

        ### If three ranges forward from i0 have positive ranges, ...
        if min(ranges[i0:i0+3]) > 0e0:
          ### ... and middle range is less than or equal to end ranges ...
          if min(ranges[i0:i0+3:2]) >= ranges[i0+1]:
            iMin = i0 + 1                         ### ... use middle range

        ### Next increment forward or backward in time
        i += 1

      """
      print(dict(i0=i0,iMin=iMin,iCenter=iCenter,ltcorrIn=ltcorrIn
                ,utcEst=sp.etcal(etEst),dt=dt
                ,utcplus=sp.etcal(etEst+dt)
                ,utcminus=sp.etcal(etEst-dt)
                ))
      """

      ### If NO success at any time and ltcorrIn was not specified
      ### (so FTCA_GETRANGE defaulted to 'LT' correction), then try
      ### repeating the calculations without the correction

      self.msg = 'OK'
      if iMin < 0:

        ### Fail with an error message if ltcorrIn argument is 'NONE' ...

        if ltcorrIn.upper() == 'NONE':
          self.msg = 'Could not find TCA'
          assert False

        ### else call self recursively with no light time correction

        tmp = FINDTCA(targ, obs, utcEstArg=utcEstArg, ltcorrIn='NONE')

        ### Copy recursed instance contents to current self instance
        self.__dict__ = tmp.__dict__

        ### Save recursed message, write new message
        self.none_ltcorr_msg = self.msg
        self.msg = 'WARNING:  FAILED WITH LIGHT TIME CORRECTION; RE-TRYING WITHOUT ...'

        #assert False,self.msg

      ### Use Newton-Raphson to find TCA

      etEst = etEst + (iMin-iCenter) * delEst
      ranje = ranges[iMin]
      delEstActual = delEst

      i=0
      while i <= 10 and delEst != 0e0 and delEstActual != 0e0:
        oldRange = ranje
        oldEt = etEst
        (ranje,stOut,unused1,unused2
        ,) = ftca_getrange(targ,obs,oldEt,ltcorrIn=self.lt,debug=True)
        p=stOut[:3]
        v=stOut[3:]
        try: delEst = - sp.vdot(p,v) / sp.vdot(v,v)
        except:
          traceback.print_exc()
          print(dict(v=v,p=p,targ=targ,obs=obs,ftca_getrange=(ranje,stOut,unused1,unused2,)))
          raise
        etEst = oldEt + delEst
        delEstActual = etEst - oldEt
        i=i+1

      if delEstActual != 0e0: oldRange = ranje

      self.et = abs(etEst)
      self.tca_err = abs(delEst)
      self.ca_dist_err:  abs(ranje-oldRange)

      (ranje,stOut,ltimeOut,self.lt_corr
      ,) = ftca_getrange(targ,obs,etEst,ltcorrIn=self.lt)

      targPosJ2k = stOut[:3]
      targVelJ2k = stOut[3:]
      self.obs_pos_j2k = sp.vscl(-1e0,targPosJ2k)
      self.obs_vel_j2k = sp.vscl(-1e0,targVelJ2k)
      self.miss_distance = sp.vnorm(stOut[:3])
      self.speed = sp.vnorm(stOut[3:5])

      ### Define trajectory plane reference frame as:
      ###   XY plane contains obs to targ and obs velocity
      ###   +X = vector to target from observer at CA
      ###   +Y near observer velocity

      mtx_trajPln = sp.twovec(targPosJ2k,1,self.obs_vel_j2k,2)

      (sunRange,stSun,unused1,unused2
      ,) = ftca_getrange('SUN',targ,etEst-ltimeOut,ltcorrIn=self.lt)

      vld = sunRange > -1e0   ### Sun lookup was valid

      self.sun_pos_j2k = stSun[:3] if vld else self.obs_vel_j2k
      self.sun_status = vld and 'OK' or 'SPICE CALL FAILED; USING UPTRACK'

      self.sun_azel_fb = sp.vscl(sp.dpr(),sp.reclat(sp.mxv(mtx_trajPln,self.sun_pos_j2k)))[1:]

      self.JDTdbTca = sp.j2000() + (etEst / spd)
      etcal = sp.etcal(etEst)
      self.TDBTIMEIDL = f'{etcal} TDB (NOT IDL CALDAT)'
      self.CALTIMESPICE = f'{etcal} TDB (SPICE ETCAL)'

      try:
        self.UTCTIME = sp.et2utc(etEst,'C',3,999) +' UTC (SPICE ET2UTC)'
        timout_fmt = 'YYYY MON DD HR:MN:SC.### TDB (SPICE TIMOUT)::TDB'
        self.TDBTIMESPICE = sp.timout(etEst,timout_fmt,999)
      except:
        self.TDBTIMESPICE = self.UTCTIME = 'LEAPSECOND-Kernel missing?'

      self.tofCalc = findvinf.FINDVINF(self
                                      ,VinfTarg
                                      ,VinfEtOffset,VinfEtStep,VinfEtNum
                                      )

    except:
      traceback.print_exc()
      self.msg = traceback.format_exc()

    ### Clean up kernels if requested
    try:
      if not nounloadkernels: ftca_furnsh(kernels,unload=True)
    except: pass


"""
########################################################################
### FTCA_TEST:  FINDTCA support function
#######################################################################n
function ftca_test, targ, obs, kernels=kernels
             , VinfTarg=VinfTarg
             , debug=debug

  ftca_furnsh, kernels, debug=debug

  TcaOut = findtca(targ, obs, VinfTarg=VinfTarg, debug=debug)

  ftca_furnsh, kernels, /unload, debug=debug

  return, TcaOut

end

bsp='nh_ref_v_od059b_encounter_only.bsp'
tls='naif0010.tls'

testEr=0
catch,testEr
if testEr == 0 then begin
  caOut=ftca_test('999',kernels=[bsp],VinfTarg='999',/debug)
  catch,/cancel
endif else begin
  catch,/cancel
  message,/continue,!error_state.msg
  message,/reset
endelse
help,/st,caOut
end
"""
