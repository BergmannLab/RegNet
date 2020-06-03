clear inopts
close all

% Function
funName = 'fkjellstrom2'

i=1;

%outKjell2_2 = cell(1,25);

% Dimension
for N = 25

    % From simulation
    stopKjell=fkjellstrom2(2.34861543*ones(N,1));
    
    % number of iterations
    for j=1:1%25
        tic
        j
        % Set-up of GaA options
        inopts.VerboseModulo = N*100;
        inopts.bSaving='on';
        inopts.SavingModulo=N*10;
        inopts.mode = 1;
        inopts.LBounds = zeros(N,1);
        inopts.UBounds = 2*pi*ones(N,1);
        inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
        inopts.MaxIter = 1e4*N;
        inopts.StopFitness = stopKjell+1e-6;
        inopts.valP = 1/exp(1);
        inopts.Display = 'off';
        inopts.Plotting = 'on';
        inopts.N_C=1e6;
        inopts.N_T=1;
        inopts.N_mu=1;
        inopts.beta=0.1;
        inopts.gRestart = 1;
        
        % Start value
        xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);

        [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
        %outKjell2_2{i,j} = out;
        toc
    end
    i=i+1;
end