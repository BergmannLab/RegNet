% README for MATLAB Gaussian Adaptation (GaA) code
%
% Christian L. Mueller
% MOSAIC group, ETH Zurich, Switzerland
%
% Prior versions of this code have been used to produce the results in:
%
% C. L. Mueller and I. F. Sbalzarini. Gaussian Adaptation revisited - an
% entropic view on Covariance Matrix Adaptation.
% In Proc. EvoStar, volume 6024 of Lecture Notes in Computer Science,
% pages 432?441, Istanbul, Turkey, April 2010. Springer.
%
% C. L. Mueller and I. F. Sbalzarini. Gaussian Adaptation as a unifying
% framework for continuous black-box optimization and adaptive Monte Carlo sampling.
% In Proc. IEEE Congress on Evolutionary Computation (CEC), Barcelona,
% Spain, July 2010.
%
% These papers can also be used for an introduction to the working principles of GaA.

I briefly describe the files contained in this package

a) Algorithm

gaussAdapt.m

This file contains the Gaussian Adaptation code. This code is in-line commented. It is inspired by
Niko Hansen's MATLAB implementation of CMA-ES und similar to use. The function reads:

% function [xmin,fmin,counteval,out] = gaussAdapt(fitfun,xstart,inopts)
% Input: 
% fitfun: Name of the fitness/target function as string or function handle
% xstart: initial candidate solution/sample
% inopts: option structure that determines internal strategy parameters (see code for details)
%
% Output:
% xmin: minimum candidate solution found by GaA (when using GaA as optimizer)
% fmin: fitness value of the xmin (when using GaA as optimizer)
% counteval: Number of function evaluations
% out: Output structure storing all relevant information (see code for details)

gaussAdapt also contains several child functions (from line 874 on). Most of them are test functions 
for optimization, design centering and sampling. 

b) Test scripts

All files starting with testGaA* contain scripts that show how to initialize and call GaA for the various scenarios.

1.) Design centering (inopts.mode = 0) : 

testGaAConRegion.m:    Shows how GaA approximates the feasible region of a linearly constrained 2D region  
testGaAConElli.m:          Shows how GaA approximates the feasible region of an nD ellipsoidal region. It also contains 
                                      details about how to compute the approximate volume of the underlying region. 

2.) Optimization (inopts.mode = 1) :

testGaASphere.m:          Sphere function
testGaARosen.m:           Rosenbrock's function
testGaAMullerBrown.m: Muller-Brown 2D landscape
testGaADfunnel.m:        Lunacek's double funnel function
testGaAKjellstrom2.m:  Kjellstrom's function 2
testGaARast.m:             Rastrigin's function
testGaANoisyS.m:         Noisy sphere function (additive noise)
testGaACEC2005.m:      If the CEC 2005 benchmark suite is available this script can be used to run GaA on it

3.) MCMC sampling (inopts.mode = 2) :

testGaASampler.m:        Shows how to sample from several target distributions (Liang's and Haario's test cases)
testGaANeal.m:             Shows how to sample from Neal's funnel distribution (a case where M-GaA should fail!)

c) Support files:

sphereVol.m:                Computes the volume of a nD sphere (used for volume computation)
error_ellipse.m:            Plots a 2D/3D Gaussian distributions (used for displaying the search trajectory)
LiangExMat.mat:           Data file for one of the target test distributions

 
IN NO EVENT SHALL THE ETH BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE ETH HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE ETH SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE ETH HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.



 


