import pbe2mdl
import spiceypy as sp
sp.furnsh('a.tm')
ucs = """
First user comment
Second user comment
...
Last user comment
""".strip().split('\n')
uls = """
First user label
Second user label
...
Last user label
""".strip().split('\n')
pbe2mdl.pbe2mdl([10,2,5],5,"Filename.ext"
               ,userCommentsArg=ucs
               ,userLabelArg=uls
               ,resolutionArg=45
               )
