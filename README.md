# CFMM-PB
This Repo publishes and maintains the source code for A Cartesian FMM-accelerated Boundary Integral Poisson-Boltzmann Solver developed by Jiahui Chen, Johannes Tausch, and Weihua Geng.

Reference: https://arxiv.org/abs/2110.13778

Acknowlegement: The work of W.G. and J.C. was supported by NSF grant DMS-1418957, DMS-1819193, SMU new faculty startup fund and SMU center for scientific computing. The work of J.T. is in part funded by NSF grant DMS-1720431.

## Usage
1. Install the triangular surface generators MSMS (https://ccsb.scripps.edu/msms/downloads/) and/or NanoShaper (https://gitlab.iit.it/SDecherchi/nanoshaper); Make sure the binary executive file of one or both were added into the system path (e.g. modify .bashrc file for Linux or .bash_profile file for MacOS)

2. Clone the code and run Makefile to build executive file coulomb

3. $ ./coulomb proteinname [options]

4. proteinname.pqr should be copy to the folder test_proteins.

5. [options] are specified as below (one might specify the values followed "=" to replace the default value as shown):

    -S=0.8
    the percentage of how much FMM is used

    -o=1e-6
    tolerance for GMRES

    -p=-1
    Taylor expansion order, negative with adaptive FMM

    -q=1
    Quadrature order

    -t=5
    tree level

    -d=1
    MSMS density

    -eps1=1 
    -eps2=80
    dielectric constant eps1 inside, eps2 outside

    -k=0.1275
    kappa (the inverse ionic strength)

    -m=1 (surface type)
    1 chooses MSMS
    2 chooses NanoShaper





