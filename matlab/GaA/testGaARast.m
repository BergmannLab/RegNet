clear inopts
close all

% Function
funName = 'frastrigin'

NVec = [2,5,10:20:50];

numIter = 10;

outRastCell=cell(length(NVec),numIter);
i=1;

% Dimension
for N = NVec
    % number of iterations
    for j=1:numIter
        disp(['Dimenison: ',num2str(N), ', Iteration: ',num2str(j)])
        % Set-up of GaA options
        inopts.VerboseModulo = 2;
        inopts.bSaving='on';
        inopts.SavingModulo=1;
        inopts.mode = 1;
        inopts.LBounds = -5*ones(N,1);
        inopts.UBounds = 5*ones(N,1);
        inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
        inopts.StopFitness = 1e-9;
        inopts.Display = 'on';
        inopts.Plotting = 'on';
        inopts.valP = 1/exp(1);
        inopts.N_mu = 1;
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);

        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
        outRastCell{i,j}=out;
        close gcf

    end
    i=i+1;
end