import numpy
import datetime
import collections
import spiceypy as sp

try: clight
except:
  clight = sp.clight()


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
      import traceback
      traceback.print_exc()


#######################################################################n
def ftca_getrange(targ, obs, et, ltcorr='LT'):
  """
ftca_getrange:  FINDTCA support function, get range from observer (obs)
                to target (targ) at et (TDB) using light time correction
                (ltcorr=).  Also return state vector (float[6]) from the
                observer to the target, and light time

"""
  ltime = -1d0 / cspice_clight()
  try:
    cspice_spkezr,targ,et,'J2000',ltcorr,obs,stOut,ltime
    return ltime*cspice_clight(),stOut,ltime
  except:
    import traceback
    traceback.print_exc()
  return -1e0,sp.vequ([0]*6,6),-1e0/clight

def ftca_getopt(defawlt, arg, poolvars):
  """
FTCA_GETOPT:  FINDTCA support function, select option
              via argument choices in decreasing priority:
              ARG       Function argument
              POOLVARS  SPICE kernel pool variable(s)
              DEFAWLT   Default value

              return value,which

"""

  if not (arg is None): return,Arg,'ARGUMENT'

  if not isinstance(poolvars,collections.Sequence):
    return ftcs_getopt(defawlt,arg,[poolvars])

  typ=isinstance(defawlt,collections.Sequence
                ) and type(default[0]) or type(defawlt)

  for poolvar in poolvars:
    try:
      if typ is str    : Vals = sp.gcpool(poolvar,0,999,999)
      elif typ is float: Vals = sp.gdpool(poolvar,0,999)
      else             : Vals = sp.gipool(poolvar,0,999)

      if len(Vals) eq 1L then Vals = Vals.pop()
      return, Vals,'KERNELPOOL'
    except:
      continue

  return,defawlt,'DEFAULT'

"""

########################################################################
########################################################################
### FINDTCA:  Main function
########################################################################
########################################################################
function findtca, targArg, obsArg, utcEst=utcEstArg
                , VinfTarg=VinfTarg
                , VinfEtOffset=VinfEtOffset
                , VinfEtStep=VinfEtStep
                , VinfEtNum=VinfEtNum
                , ltcorrIn=ltcorrArg
                , kernels=kernels
                , nounloadkernels=nounloadkernels
                , debug=debug
  dodebug = keyword_set(debug)

  ftca_furnsh,kernels,debug=debug

  targ = ftca_getopt('PLUTO', targArg, 'FINDTCA_TARGET')
  obs = ftca_getopt('-98', obsArg, 'FINDTCA_OBSERVER',which=wObs)

  ### use UTC estimate if given as argument,
  ### else use kernel pool variable FINDTCA_UTCEST
  ### else use NOW

  defawlt = strcompress( 'JD ' + string(julday(),f='(f30.1)') ) + ' TDB'
  utcEst = ftca_getopt(defawlt, utcEstArg, 'FINDTCA_UTCEST',which=wUtcEst)

  if wUtcEst eq 'DEFAULT' then begin
    delEst=7d0
  endif else begin
    delEst=1d0
  endelse

  ### Find out if LEAPSECONDS have been loaded

  lskEr=0L
  catch,lskEr
  if lskEr eq 0L then begin
    ###cspice_deltet,0d0,'ET',d
    ###catch,/cancel
    cspice_str2et,utcEst,etEst
    catch,/cancel
  endif else begin
    catch,/cancel
    cspice_tparse,utcEst,etEst ,tpErrmsg
    message,/continue,'WARNING:  NO LEAPSECOND KERNEL AVAILABLE'
    if strtrim(tpErrmsg,2) ne '' then message,/continue,'TPARSE:  ' + tpErrmsg
  endelse

  ### Convert UTC estimate to ET (s past J2000), delEst to s

  delEst = delEst * cspice_spd()

  ### Go 250 steps either side of estimate, take first local minimum range

  iCtr=1000L
  nRanges = iCtr * 2L + 1L
  ranges = make_array(nRanges,val=-1d0)

  iMin = -1L
  i = 0L
  while i le iCtr and iMin lt 0L do begin
    i0 = iCtr+i
    ranges[i0] = ftca_getrange(targ,obs,etEst+(delEst*i),ltcorrIn=ltcorrArg)
    if min(ranges[i0-2L:i0]) gt -1d0 then begin
      if min(ranges[[i0-2L,i0]]) ge ranges[i0-1L] then iMin = i0-1L
    endif
    if iMin lt 0L then begin
      i0 = iCtr-i
      ranges[iCtr-i] = ftca_getrange(targ,obs,etEst-(delEst*i)
                                    ,ltcorrIn=ltcorrArg)
      if min(ranges[i0:i0+2L]) gt -1d0 then begin
        if min(ranges[[i0+2L,i0]]) ge ranges[i0+1L] then iMin = i0+1L
      endif
    endif

    i = i + 1L
  endwhile

  ### If NO success at any time and ltcorrArg was not specified
  ### (so FTCA_GETRANGE defaulted to 'LT' correction), then try
  ### repeating the calculations without the correction

  if iMin lt 0L then begin

    ### Fail with an error message if ltcorrArg was specified ...

    if n_elements(ltcorrArg) eq 1L then message,'Could not find TCA'

    ### else call self recursively with no light time correction ('NONE')

    message,/continue
           ,'WARNING:  FAILED WITH LIGHT TIME CORRECTION# RE-TRYING WITHOUT ...'

    rtn = findtca( targ, obs, utcEst=utcEstArg, ltcorrIn='NONE', debug=debug)
    if not keyword_set(nounloadkernels) then ftca_furnsh,kernels,/unload
    return, rtn

  endif

  ### Use Newton-Raphson to find TCA

  etEst = etEst + (iMin-iCtr) * delEst
  range = ranges[iMin]
  md = 0d0
  delEstActual = delEst

  i=0L
  while i le 10L and delEst ne 0d0 and delEstActual ne 0d0 do begin
    oldRange = range
    oldEt = etEst
    range = ftca_getrange(targ,obs,oldEt,stOut=stOut,ltcorrIn=ltcorrArg)
    p=stOUt[0:2]
    v=stOut[3:5]
    delEst = - cspice_vdot(p,v) / total(v*v)
    etEst = oldEt + delEst
    delEstActual = etEst - oldEt
    if dodebug then begin
      oldMd = md
      cspice_vperp, p, v, mdVec
      md = cspice_vnorm(mdVec)
      if i eq 0L then begin
        print, 'Iter', 'deltaEtCalc', 'deltaEtActual'
             , 'deltaMissDist', 'MissDist'
             , f='(a4,4a16)'
      endif else begin
        print,i,delEst,etEst-oldEt,md-oldMd, md
             , f='(i4,4g16.8)'
      endelse
    endif
    i=i+1
  endwhile

  if delEstActual ne 0d0 then oldRange = range
  range = ftca_getrange(targ,obs,etEst,stOut=stOut,ltcorrIn=ltcorrArg
                       ,ltcorrOut=ltcorrOut
                       ,ltimeOut=ltimeOut)

  targPosJ2k = stOut[0:2]
  targVelJ2k = stOut[3:5]
  obsPosJ2k = -targPosJ2k
  obsVelJ2k = -targVelJ2k

  ### Define reference frame as:
  ###   XY plane contains obs to targ and obs velocity
  ###   +X = vector to target from observer at CA
  ###   +Y near observer velocity

  cspice_ucrss, targPosJ2k, obsVelJ2k, trajPlaneZ     ### +X cross ~+Y = +Z
  cspice_ucrss, trajPlaneZ, targPosJ2k, trajPlaneY    ### +Z cross +X  = +Y
  cspice_ucrss, trajPlaneY, trajPlaneZ, trajPlaneX    ### +Y cross +Z  = +X

  sunRange = ftca_getrange('SUN',targ,etEst-ltimeOut,stOut=stSun
                       ,ltcorrIn=ltcorrArg
                       ,ltimeOut=ltimeOut)

  if sunRange gt -1d0 then begin
    sunPosJ2k = stSun[0:2]
    sunStatus = 'OK'
  endif else begin
    sunPosJ2k = obsVelJ2k
    sunStatus = 'SPICE CALL FAILED# USING UPTRACK'
  endelse

  sunAzElRFb = cv_coord( /degree, /to_sph
                       , from_rec=[ cspice_vdot(trajPlaneX,sunPosJ2k)
                                  , cspice_vdot(trajPlaneY,sunPosJ2k)
                                  , cspice_vdot(trajPlaneZ,sunPosJ2k)
                                  ]
                       )

  jdTdbTca=cspice_j2000()+etEst/cspice_spd()
  caldat,jdTdbTca,month,day,year,hour,minute,seconds
  month=strmid('JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC',month*3-3,3)
  TdbTimeIdl = string(year,' ',month,' ',day,f='(i4.4,a,a,a,i2.2)' )
               + ' '
               + string( hour, ':',minute,':',seconds, f='(2(i2.2,a),f06.3)')
               + ' TDB (IDL CALDAT)'

  cspice_etcal,etEst,CalTimeSpice

  utcEr=0L
  catch,utcEr
  if utcEr eq 0L then begin
    cspice_et2utc,etEst,'C',3,UtcTimeBare
    cspice_timout,etEst,'YYYY MON DD HR:MN:SC.### TDB (SPICE TIMOUT)::TDB',999
                ,TdbTimeSpice
    catch,/cancel
    UtcTime=UtcTimeBare + ' UTC (SPICE ET2UTC)'
  endif else begin
    catch,/cancel
    UtcTime = 'LEAPSECOND-Kernel missing?'
    TdbTimeSpice = UtcTime
  endelse

  xxx =
  {            ET:  etEst
  , MISS_DISTANCE:  cspice_vnorm(stOut[0:2])
  ,         SPEED:  cspice_vnorm(stOut[3:5])
  ,        TARGET:  targ
  ,           OBS:  obs
  ,   OBS_POS_J2k:  obsPosJ2k
  ,   OBS_VEL_J2k:  obsVelJ2k
  ,   SUN_POS_J2k:  sunPosJ2k
  ,   SUN_AZEL_FB:  sunAzElRFb[0:1]
  ,    SUN_STATUS:  sunStatus
  ,       LT_CORR:  ltcorrOut
  ,       TCA_ERR:  abs(delEst)
  ,   CA_DIST_ERR:  abs(range-oldRange)
  ,      JDTDBTCA:  jdTdbTca
  ,    TDBTIMEIDL:  TdbTimeIdl
  ,  CALTIMESPICE:  CalTimeSpice + ' TDB (SPICE ETCAL)'
  ,  TDBTIMESPICE:  TdbTimeSpice
  ,       UTCTIME:  UtcTime
  }

  tcaOut = create_struct( xxx
  , 'tofCalc', findvinf( xxx
                       , VinfTarg
                       , VinfEtOffset
                       , VinfEtStep
                       , VinfEtNum
                       , debug=debug
                       )
  )

  if not keyword_set(nounloadkernels) then ftca_furnsh,kernels,/unload

  return, tcaOut
end
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

testEr=0L
catch,testEr
if testEr eq 0L then begin
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
