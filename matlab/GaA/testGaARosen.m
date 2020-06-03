clear inopts
close all

% Function
funName = 'frosen'

NVec = [10,20,30];

numIter = 1;

outRosenCell=cell(length(NVec),numIter);
i=6;

% Dimension
for N = NVec
    
    % number of iterations
    for j=1:numIter
        
        disp(['Dimenison: ',num2str(N), ', Iteration: ',num2str(j)])
        
        inopts.mode = 1;
        
        % Set-up of GaA options
        inopts.VerboseModulo = 1000;
        inopts.bSaving='on';
        inopts.SavingModulo=10;
        
        inopts.LBounds = -2.048*ones(N,1);
        inopts.UBounds = 2.048*ones(N,1);
        inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
        inopts.StopFitness = 1e-9;
        
        inopts.MaxIter=N*1e4;
        inopts.Display = 'off';
        inopts.Plotting = 'on';
        inopts.valP = 1/exp(1);
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);

        outRosenCell{i,j}=out;
        close gcf
    end
    i=i+1;
end

