import os
import pprint
import pbecalcs

os.system("kernels/meta_kernel.tm --get-kernels")

pc = pbecalcs.PBESTRUCT(kern=["kernels/meta_kernel.tm"])
pprint.pprint((vars(pc)
              ,vars(pc.kinetx)
              ,vars(pc.pbeArr[0,0])
              ,vars(pc.pbeArr[0,0].scanVinfUptrack)
              ,))
