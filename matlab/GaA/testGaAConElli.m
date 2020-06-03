clear inopts
close all

% Function
funName = 'fConElli'

% (x-b)'A(x-b)<1 defines the feasible region, 
% i.e. an ellipsoid centered around b

% Dimensions
NVec = 2:2:20;

% Help variables
lastR = zeros(length(NVec),1);
trueVol=zeros(length(NVec),1);
approxVol=zeros(length(NVec),1);

% Counter
i=1;

% Dimension
for N = NVec
        
    disp(['Dimension: ',num2str(N)])

    % A defines an ellipsoid with decreasing semi-axis
    % A = diag(1:N);
    
    % A defines a scaled sphere
    A = 1/2.*eye(N,N);
    
    % b is the shift
    b=0.5*ones(N,1);
    
    funArgs = [A,b];
    
    % Set-up of GaA options
    inopts.VerboseModulo = 1e4;
    inopts.bSaving='on';
    inopts.SavingModulo=5e3;
    inopts.bSaveCov = 0;
    % Design centering
    inopts.mode = 0;
    inopts.LBounds = -4*ones(N,1);
    inopts.UBounds = 4*ones(N,1);
    inopts.r = 0.5;
    inopts.N_m = N*1e1;
    inopts.N_C = N*1e2;
    inopts.Display = 'off';
    inopts.Plotting = 'on';
    inopts.valP = 1/exp(1);
    inopts.c_T = 0.5;
    inopts.funArgs = funArgs;
    inopts.MaxIter = N*1e4;
    
    % Start value, must be a feasible point
    xstart = 0.55*ones(N,1);
    
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
    
    % True volume of the ellipsoid
    trueVol(i) = det(inv(sqrt(A)))*sphereVol(1,N);
    
    % Extract volume information
    if strcmp(inopts.bSaving,'on')
        lastR(i) = mean(out.rVec(end));
        if inopts.bSaveCov == 1
            C = lastR(i)^2*out.Q{end}'*out.Q{end};
        else
            C = lastR(i)^2*out.lastQ'*out.lastQ;
        end
    else
        lastR(i) = out.lastR;
        C = lastR(i)^2*out.lastQ'*out.lastQ;
    end
    
    % Just for checking
    % detC=det(C) % should be == lastR(i)^(2*N)
    
    % Scaling according to chi2 distribution
    c_p=chi2inv(out.P_emp,N);
    
    % Do an eigendecomposition of the true covariance
    [eigVecs,eigVals] = eig(C);
    
    % Rescale the covariance
    elliAxes =  sqrt(c_p*eigVals);
    
    % Compute the approximate volume of the minimum volume ellipsoid that
    % corresponds to the statistical description
    approxVol(i)=prod(diag(elliAxes))*sphereVol(1,N);
    
    trueVol(i)
    approxVol(i)
    i=i+1;
    
end

figure
semilogy(NVec,trueVol)
hold on
semilogy(NVec,approxVol,'r')
xlabel('dimension')
ylabel('Volume')
grid on
title('True volume(blue) and approximate volume (red)')











