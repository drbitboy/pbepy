import spiceypy as sp
from kinetx_mu69_20180626 import KINETX_MU69_20180626
########################################################################
### Wrapper for current info from Kinetx
########################################################################
### 20180626:  Removed Pluto
def kinetxcurrent(target=None):

  sTarget = str(2486958 if target is None else target).strip()

  try: iTarget = sp.bods2c(sTarget)
  except:
    import traceback
    traceback.print_exc()
    assert False,f'Could not understand target [{sTarget}]'

  ### MU69 => SPICE ID 2486958
  ### MU69_BARYCENTER => SPICE ID 2486959
  ### - Multiple bodies SPICE IDs are 218695801, 148695802, etc.
  if iTarget == 2486958 or (iTarget//100) == 2486958:
    return KINETX_MU69_20180626()    ### 2018-06-26

  assert False, f'Target {sTarget} is unrelated to MU69'
