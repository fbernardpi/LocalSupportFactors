***************************************************************************

Overview:

This is code to do structured matrix factorizations based on the
formulation given in the "Applications" section of the paper:

Haeffele, B., Young, E., and Vidal, R. "Structured low-rank matrix 
factorization: Optimality, algorithm, and applications to image 
processing." ICML. 2014.

The code here is a slightly generalization of the formulation given in the
paper as it also includes options to include additional 'non-factorized'
variables in the optimization (for example, to add intercepts).
***************************************************************************

Files:

FactorTVL1L2_v1.m - The main function.

ToyExample.m - A script with simulated data to show how to run the code.

proximalTVL1L2.m - Function to calculate the proximal operator of the sum
    of the L1, L2, and Total-Variation norms, as given in equation (19) of
    the above paper.

mexProximalL1TV.cpp - The proximalTVL1L2.m function calls this mex file to 
    calculate the proximal operator of the L1+TV term.  The calculation is
    done by doing coordinate descent on the dual problem and checking the
    duality gap for the stopping criteria.  See below for instuctions on
    compiling.

L1TVprox.cpp/h - Library and header files for the mexProximalL1TV.cpp mex
    function.

GetLatticeIndex.m - Returns the pairs of neighboring indices for an
    8-connected 2D lattice (i.e. neighboring pixels of an image of a 
    given size).

MakeFourierD.m - Makes a dictionary of Fourier components as an example for
    the ToyExample.m script.
***************************************************************************

Compiling the mex file:

The implimentation was purposely written for ease of compilation (i.e. no 
external libraries, such as LAPACK, are required).  Assuming the mex C/C++
compiler is properly configured and that all of the source files are in the
current directory, the mex file should compile with the following command:

mex mexProximalL1TV.cpp L1TVprox.cpp
***************************************************************************

Contact:

Please send any questions or found bugs to Ben Haeffele
bhaeffele@jhu.edu