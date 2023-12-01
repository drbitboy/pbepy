#!/usr/bin/env python

"""
\begindata

PATH_VALUES  = ( './kernels' )
PATH_SYMBOLS = ( 'KERNEL_DIR' )

KERNELS_TO_LOAD = (
  '$KERNEL_DIR/naif0012.tls'
  '$KERNEL_DIR/nh_recon_arrokoth_od147_v01.bsp'
  '$KERNEL_DIR/nh_pcnh_008.tpc'
  '$KERNEL_DIR/pck00010.tpc'
)
\begintext

Notes

For testing Git repository pbepy

Notes
=====

* Edit PATH_VALUES if kernels are elsewhere.
* Retrieve the kernels, run this script:

    % ./kernels/meta_kernel.tm --get-kernels

"""

import os
import sys

if "__main__" == __name__ and ['--get-kernels'] == sys.argv[1:2]:

  equpar = '= ('.split()
  plusequpar = '+= ('.split()
  mkprefix = '$KERNEL_DIR/'
  Lmkprefix = len(mkprefix)
  looking_for_path = True
  looking_for_kernels = False

  kernel_list = list()

  local_path = None

  for docline in [s.strip()
                  for s in __doc__.replace("'","").split('\n')
                  if s.strip()
                 ]:

      toks = docline.split()
      if not toks: continue
      tok0 = toks.pop(0)
      tokcount = len(toks)

      if looking_for_kernels:

        if not tokcount:

          ### $KERNEL_DIR/base.ext
          if tok0.startswith(mkprefix):
            kernel_list.append(tok0[Lmkprefix:])

          ### End of KERNELS_TO_LOAD
          if ')' == tok0: looking_for_kernels = False

        continue

      if 'KERNELS_TO_LOAD' == tok0 and 2 == tokcount:
        toks2 = toks[:2]
        if equpar == toks2 or plusequpar == toks2:
          looking_for_kernels = True
          continue

      if looking_for_path and 4 == tokcount:

        if 'PATH_VALUES' == tok0 and ')' == toks[-1]:
          toks2 = toks[:2]
          if equpar == toks2 or plusequpar == toks2:
            local_path = toks[2]
            looking_for_path = False
            continue

  wprefix = 'https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data'
  ext2subdir = lambda k: k.split('.').pop()[1:] + 'k'

  for kernel in kernel_list:
    if os.path.exists(os.path.join(local_path,kernel)):
      print(f'Skipping {kernel} ...')
      continue

    if not os.path.exists(local_path): os.system(f'mkdir -pv {local_path}')

    wpath = os.path.join(wprefix, ext2subdir(kernel), kernel)
    print(f'Retrieving [{wpath}] ...')
    os.system(f'wget -q {wpath} -P {local_path}')

  print( '... Done')
"""
wget https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/lsk/naif0012.tls -P ./kernels
wget https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/nh_recon_arrokoth_od147_v01.bsp -P ./kernels
wget https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/pck/nh_pcnh_008.tpc -P ./kernels
wget https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/pck/pck00010.tpc -P ./kernels
"""
