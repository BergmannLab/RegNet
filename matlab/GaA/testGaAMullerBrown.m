clear inopts
close all

% Function
funName = 'MullerBrown'

% Dimension
N = 2

% From simulation
stopMuller= -151.666;

% number of iterations
for j=1:1%25
    j
    % Set-up of GaA options
    inopts.VerboseModulo = N*100;
    inopts.bSaving='on';
    inopts.SavingModulo=1;
    inopts.mode = 1;
    inopts.LBounds = [-1.5;-0.9]
    inopts.UBounds = [1.1;2.1]
    inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
    inopts.MaxIter = 1e4*N;
    inopts.StopFitness = stopMuller+1e-6;
    inopts.valP = 1/exp(1);
    inopts.Display = 'on';
    inopts.Plotting = 'on';
    
    % Start value
    xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
    
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
end
