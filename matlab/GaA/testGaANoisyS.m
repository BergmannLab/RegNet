clear inopts
close all

% Function
funName = 'fnoisysphere'

stdVals = [1:4];
numIter=1;

outAddNoisyS = cell(length(stdVals),numIter);

% Dimension
N = 30

i=1;
for s=stdVals
    
    close gcf

    s
    % number of iterations
    for j=1:numIter
        tic
        j
        % Set-up of GaA options
        inopts.mode = 1;
        
        inopts.VerboseModulo = 1e2*N;
        inopts.bSaving='on';
        inopts.SavingModulo = 1e1*N;
        inopts.LBounds = -5*ones(N,1);
        inopts.UBounds = 5*ones(N,1);
        inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
        inopts.MaxIter = 3e3*N;
        inopts.StopFitness = -Inf;
        inopts.valP = 1/exp(1);
        
        inopts.Display = 'on';
        inopts.Plotting = 'on';
        inopts.bRestart = 0;
        
        inopts.TolR = 1e-20;
        inopts.funArgs = s;
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
        
        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
        outAddNoisyS{i,j} = out;
        toc
    end
    
    i=i+1;
end