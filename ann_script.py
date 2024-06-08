#!/usr/bin/env python3

"""Import PBE modules"""
import pbecalcs
import pbe2mdl

"""Instantiate  PBE structure for default Arrokoth/NH flyby"""
pbec=pbecalcs.PBESTRUCT(kern=['kernels/meta_kernel.tm'],sigm=[100,66,33],nSig=3)

"""Write MDL file from that PBE structure"""
pbe2mdl.pbe2mdl(pbec,1.0,'testing/ann240607.mdl',mtxJ2k2UncertArg=pbec.mtx_j2k2Uncert)
