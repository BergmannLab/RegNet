% Script that runs Gaussian Adaptation on the sphere function
clear inopts
close all

% Function
funName = 'fsphere'

NVec = [2,5,10:10:50];

numIter = 10;

outSphCell=cell(length(NVec),numIter);
i=1;

% Dimension
for N = NVec
    
    % number of iterations
    for j=1:numIter
        disp(['Dimenison: ',num2str(N), ', Iteration: ',num2str(j)])
        
        % Set-up of GaA options
        inopts.mode = 1;
        
        inopts.VerboseModulo = Inf;
        inopts.bSaving='on';
        inopts.SavingModulo=100;
        inopts.LBounds = -5*ones(N,1);
        inopts.UBounds = 5*ones(N,1);
        inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
        
        inopts.MaxIter = 1e4*N;
        inopts.StopFitness = 1e-9;
        
        inopts.Display = 'off';
        inopts.Plotting = 'off';
        inopts.valP = 1/exp(1);
        inopts.bRestart=0;
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
        
        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
        outSphCell{i,j}=out;
    end
    i=i+1;
end










