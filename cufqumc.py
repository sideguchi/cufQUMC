#!/usr/bin/env python
"""QU-fitting using replica exchange MCMC (parallel tempering)."""

"""
MIT License

Copyright (c) 2023 Shinsuke Ideguchi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import numpy as np
import os
import argparse
import ast
import sys
import glob
import subprocess
from scipy.io import FortranFile
from collections import OrderedDict


#
# Run Fortran code
#
def wrap_Fort_mcmc(DirList,Nbeta,Nburnin,Nsample,Nadjust,par,dirname,Nexchevry,Nexch,Npar,Nthreads):

    Ndata = len(DirList)

    # MCMC parameters
    rnd = np.random.RandomState(None).randint(1,2**31)
    outarray = np.array([Nbeta,Nburnin,Nsample,Nadjust,Nexchevry,Nexch,Npar,Nthreads,rnd])
    np.savetxt(dirname+'mcmc.txt',outarray,fmt='%i')

    # par
    for i in range(Ndata):
        outarray = par[:,:,:,i].reshape(4,Npar,Nbeta).T
        FortranFile(dirname+DirList[i]+'par.dat','w').write_record(outarray)

    ### run MCMC
    LibPath = os.path.dirname(os.path.abspath(__file__))
    cmd = LibPath+'/libcufqumc {}'.format(dirname)
    subprocess.call(cmd.split())


#
# model
#
par_dict = {
'init':{'RM_radm2':0.,'fracPol':0.01,'psi0_rad':0.0,'sigma_RM':0.01,'delta_RM':0.01,'delta':0.0,'pav':0.0,'pli':0.0},
'step':{'RM_radm2':10.,'fracPol':0.1,'psi0_rad':0.2,'sigma_RM':1.,'delta_RM':1.,'delta':0.1,'pav':1.0,'pli':0.1},
'min':{'RM_radm2':-2000.,'fracPol':0.,'psi0_rad':-np.pi/2.,'sigma_RM':0.,'delta_RM':0.,'delta':-1.,'pav':-10.,'pli':-10.},
'max':{'RM_radm2':2000.,'fracPol':10.,'psi0_rad':np.pi/2.,'sigma_RM':200,'delta_RM':200,'delta':1.,'pav':10.,'pli':10.}
}

par_kind = OrderedDict(
    delta=('RM_radm2','fracPol','psi0_rad'),
    delta_pli=('RM_radm2','fracPol','psi0_rad','pli'),
    gauss=('RM_radm2','fracPol','psi0_rad','sigma_RM'),
    gauss_pav=('RM_radm2','fracPol','psi0_rad','sigma_RM','pav'),
    gauss_pli=('RM_radm2','fracPol','psi0_rad','sigma_RM','pli'),
    gauss_pav_pli=('RM_radm2','fracPol','psi0_rad','sigma_RM','pav','pli'),
    tophat=('RM_radm2','fracPol','psi0_rad','delta_RM'),
    tophat_pav=('RM_radm2','fracPol','psi0_rad','delta_RM','pav'),
    tophat_pli=('RM_radm2','fracPol','psi0_rad','delta_RM','pli'),
    tophat_pav_pli=('RM_radm2','fracPol','psi0_rad','delta_RM','pav','pli'),
    osul17=('RM_radm2','fracPol','psi0_rad','sigma_RM','delta_RM'),
    skew=('RM_radm2','fracPol','psi0_rad','sigma_RM','delta'),
    skew_pli=('RM_radm2','fracPol','psi0_rad','sigma_RM','delta','pli')
)

allmodel = [
'delta',
'delta_pli',
'gauss',
'gauss_pav',
'gauss_pli',
'gauss_pav_pli',
'tophat',
'tophat_pav',
'tophat_pli',
'tophat_pav_pli',
'osul17',
'skew',
'skew_pli'
]

comp_par_num = np.array([3,4,4,5,5,6,4,5,5,6,5,5,6],dtype=np.int32)


### beta distribution
def dist_beta(Nbeta,betamin):
    if Nbeta < 1:
        print('{}: {}: error: Nbeta should be >= 1: invalid value: "{}"'.format(sys.argv[0],sys._getframe().f_code.co_name),Nbeta)
        sys.exit()
    beta = np.empty((Nbeta),dtype=np.float64)
    if Nbeta == 1:
        beta[0] = 1.
    elif Nbeta == 2:
        beta[0] = 0.
        beta[1] = 1.
    else:
        beta[0] = 0.
        for i in range (1,Nbeta):
            beta[i] = betamin**((Nbeta-2-(i-1))/(Nbeta-2))
    return beta


def main():

    # Parse the command line options
    parser = argparse.ArgumentParser(description='CUDA Fortran QUMC.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-d',dest='DataDir',type=str, default=None,
                        help='Directory containing data files.')
    parser.add_argument('-m', dest='FitModel', type=str, default=None,
                        help='Fitting Model. ["delta", "gauss", "tophat", "osul17"] are available.\ne.g. One "delta" and two "gauss" components: -m="{\'delta\':1,\'gauss\':2}".')
    parser.add_argument('-b', dest='Nburnin', type=int, required=True,
                        help='Number of chain for Burn-in (step widths adjustment) [Required].')
    parser.add_argument('-s', dest='Nsample', type=int, required=True,
                        help='Number of chain for sampling [Required].')
    parser.add_argument('-o', dest='OutputDir', type=str, default='./out',
                        help='Output directory ["./out"].')
    parser.add_argument('-n', dest='Nbeta', type=int, default=1,
                        help='The number of Replica [100]. beta=[1] when Nbeta=1.')
    parser.add_argument('-r', dest='resume', action='store_true',
                        help='Resume MCMC from the last results in "OutputDir" directory [False]. \nWhen activated, only "Nburnin", "Nsample", "OutputDir" and "Nthreads" are needed to specify.')
    parser.add_argument('-t', dest='Nthreads', type=int, default=1024,
                        help='Number of threads per block [1024].')
    args = parser.parse_args()


    ### add '/'
    OutputDir = args.OutputDir
    if OutputDir[-1] != '/':
        OutputDir += '/'

    resume = args.resume

    if not resume:
        ### add '/'
        DataDir   = args.DataDir
        if DataDir[-1] != '/':
            DataDir += '/'

        ### Fit Model
        FitModel  = args.FitModel

        ### Number of beta
        Nbeta = int(args.Nbeta)

        ### DataList
        DataList = sorted(glob.glob(DataDir+'*'))

        ### output-directory name list
        filename_list = []
        for i in DataList:
            filename_list.append(os.path.split(os.path.splitext(i)[0])[-1]+'/')
    else:
        Nbeta = int(np.loadtxt(OutputDir+'mcmc.txt')[0])


    ### try to create output directory
    if not resume:
        # outputdir
        os.makedirs(OutputDir)
        # datadir
        for i in filename_list:
            os.makedirs(OutputDir+i)
        # output datadir
        with open(OutputDir+'datalist.txt', 'w') as f:
            for i in filename_list:
                f.write('"{}"\n'.format(i))
    else:
        with open(OutputDir+'datalist.txt', 'r') as f:
            filename_list = f.readlines()
        for i,j in enumerate(filename_list):
            filename_list[i] = j.rstrip('\n')

    ### Ndata
    Ndata = len(filename_list)

    ### Model
    if not resume:
        model = ast.literal_eval(FitModel)
        src_list = np.zeros((len(allmodel)),dtype=np.int32)
        for i,j in enumerate(allmodel):
            val = model.get(j)
            if val is not None and val > 0:
                src_list[i] = val
        # output
        np.savetxt(OutputDir+'src_list.txt',src_list,fmt='%i')
    else:
        src_list = np.loadtxt(OutputDir+'src_list.txt',dtype=np.int32)

    ### parameter
    Npar = np.sum(src_list*comp_par_num)
    par = np.empty((4,Npar,Nbeta,Ndata),dtype=np.float64)
    if not resume:
        par0 = np.empty([4,Npar,Nbeta],dtype=np.float64)
        peri = np.zeros([Npar],np.int32)
        adjnum = 0
        for modnum,modname in enumerate(allmodel):
            Nsrc = model.get(modname)
            if Nsrc is not None:
                for srcnum in range(Nsrc):
                    peri[adjnum+2] = 1
                    for typenum,typename in enumerate(par_dict):
                        for parnum,parname in enumerate(par_kind[modname]):
                            for ibeta in range(Nbeta):
                                par0[typenum,parnum+adjnum,ibeta] = par_dict[typename][parname]
                    adjnum += comp_par_num[modnum]
        for i in range(Ndata):
            par[:,:,:,i] = par0
        # output peri
        np.savetxt(OutputDir+'peri.txt',peri,fmt='%i')
    else:
        for i in range(Ndata):
            filename = OutputDir+filename_list[i]+'par.dat'
            par[:,:,:,i] = FortranFile(filename).read_record(np.float64).reshape(Nbeta,Npar,4).T

    ### read data
    if not resume:
        for i in range(Ndata):
            # data
            (freq,q,u,qerr,uerr) = np.loadtxt(DataList[i]).T
            loc  = ~np.isnan(q)
            loc *= ~np.isnan(u)
            nanchan = np.ones((len(loc)),dtype=np.int_)
            nanchan[loc] = 0
            np.savetxt(OutputDir+filename_list[i]+'data.txt',np.array([q,u,qerr,uerr,nanchan]).T)
            # Define MinBeta
            # Find the maximum S/N
            tmp1 = max(np.nanmax(np.abs(q/qerr)),np.nanmax(np.abs(u/uerr)))
            # Find the maximum difference in flux
            tmp2 = max(np.nanmax(q)-np.nanmin(q),np.nanmax(u)-np.nanmin(u))
            MinBeta = np.power(tmp1*tmp1*tmp2,-2)
            mybeta = dist_beta(Nbeta,MinBeta)
            np.savetxt(OutputDir+filename_list[i]+'beta.txt',mybeta)
        # output freq
        np.savetxt(OutputDir+'freq.txt',freq)

    ### MCMC
    Nexchevry,Nexch,Nadjust = 50,50,100
    wrap_Fort_mcmc(filename_list,Nbeta,args.Nburnin,args.Nsample,Nadjust,par,OutputDir,Nexchevry,Nexch,Npar,args.Nthreads)


if __name__ == '__main__':
    main()
