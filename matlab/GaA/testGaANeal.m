% Script that runs Gaussian Adaptation on Neal's funnel distribution
clear inopts
close all

% Function from Neal's Slice sampling paper
funName = 'fNealsFunnel';

i=1;

% Dimension
N = 10;

numIter = 1;
MaxSample = 1.e5;

tempSum=0;
% number of iterations
for j=1:numIter
    j
    % Set-up of GaA options
    inopts.VerboseModulo = 1e3;
    inopts.bSaving='on';
    inopts.SavingModulo=1;
    inopts.mode = 2;
    inopts.r = 2.34^2/N;
    inopts.N_mu = 1;
    inopts.MaxIter = MaxSample;
    inopts.Display = 'off';
    inopts.Plotting = 'on';
    inopts.valP = 0.234;
    inopts.N_C=1e4;
    
    inopts.MaxCond = N;
    inopts.MaxR = 1e2*2.34^2/N;
    
    % Start value (as given by Neal)
    xstart = [0;ones(N-1,1)]
    
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
    NealsData=out.xRaw;
end







