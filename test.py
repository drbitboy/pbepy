import os
import sys
import pprint
import pbecalcs
import pbe2mdl


if "__main__" == __name__:
  os.system("kernels/meta_kernel.tm --get-kernels")

  pc = pbecalcs.PBESTRUCT(kern=["kernels/meta_kernel.tm"])
  if '--print-pbecalcs' in sys.argv[1:]:
    pprint.pprint(dict(pbecalcs=vars(pc)
                      ,pbecalcs_kinetx=vars(pc.kinetx)
                      ,pbecalcs_pbeArr00=vars(pc.pbeArr[0,0])
                      ,yw_scanVinfUptrack=vars(pc.pbeArr[0,0].scanVinfUptrack)
                      ,yx_scanVinfDowntrack=vars(pc.pbeArr[0,0].scanVinfDowntrack)
                      ,yy_scanVinfOvertrack=vars(pc.pbeArr[0,0].scanVinfOvertrack)
                      ,yz_scanVinfBelowtrack=vars(pc.pbeArr[0,0].scanVinfBelowtrack)
                      ,zw_scanVtargUptrack=vars(pc.pbeArr[0,0].scanVtargUptrack)
                      ,zx_scanVtargDowntrack=vars(pc.pbeArr[0,0].scanVtargDowntrack)
                      ,zy_scanVtargOvertrack=vars(pc.pbeArr[0,0].scanVtargOvertrack)
                      ,zz_scanVtargBelowtrack=vars(pc.pbeArr[0,0].scanVtargBelowtrack)
                      ,))

  mdlfn = (['']+[arg for arg in sys.argv[1:] if arg.startswith('--mdl-fn=')]).pop()[9:]
  if mdlfn:
    pbe2mdl.pbe2mdl(pc,1.0,mdlfn,resolutionArg=30)
    print(f'Wrote MDL file {mdlfn}')
