"""
Usage:

  python list_targets.py .../mk/*.tm .../spk/*_[19]20000617.bsp --check

  N.B. ...20000617.bsp FURNSHes Patroclus and Menoetius ephemerides

"""
import os
import sys
import kinetxcurrent
import spiceypy as sp

def go(kernels):
  list(map(sp.furnsh,[k for k in kernels if not k.startswith('--')]))
  names = sp.gcpool('NAIF_BODY_NAME',0,999,999)
  d = dict()
  for name in names:
    spid = sp.bods2c(name)
    if spid < 0: continue
    family = spid % 100000000
    if not (family in d): d[family] = dict()
    if not (spid in d[family]): d[family][spid] = set()
    d[family][spid].add(name)

  return d

def check(d):
  for family in d:
    names = d[family]
    for spid in names:
      for name in names[spid]:
        try:
          import pprint
          print(f"Success:  {name} at {kinetxcurrent.kinetxcurrent(name).Tca.UTCTIME}")
        except:
          print(f"Failure:  {name}")
        

if "__main__" == __name__:
  d = go(sys.argv[1:])

  import pprint
  pprint.pprint(d)

  if '--check' in sys.argv: check(d)
