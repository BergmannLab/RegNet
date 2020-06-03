clear inopts
close all

% Function
funName = 'fdrastrigin'

s = 0.1;
d = 2.5^2 - s*2.5^2;

i=1;

outDfunnel = cell(numIter,1);

% Dimension
for N = 10
    
    stopDfunnel = 1e-2;
    
    % number of iterations
    for j=1:numIter
        tic
        j
        % Set-up of GaA options
        inopts.mode = 1;
        
        inopts.VerboseModulo = N*200;
        inopts.bSaving='on';
        inopts.SavingModulo=N*2;
        inopts.LBounds = -5*ones(N,1);
        inopts.UBounds = 5*ones(N,1);
        inopts.BoundPenalty = 1;
        inopts.r = 1/exp(1)*(max(inopts.UBounds)-min(inopts.LBounds));
        inopts.MaxIter = 1e4*N;
        inopts.StopFitness = stopDfunnel;
        inopts.valP = 1/exp(1);
        inopts.Display = 'on';
        inopts.Plotting = 'on';
        inopts.funArgs = [s;d];
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
        
        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
        outDfunnel{j} = out;
        toc
    end
    i=i+1;
end


