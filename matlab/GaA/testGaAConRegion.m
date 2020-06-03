clear inopts
close all

% Function
funName = 'fConRegion'

% Dimension
N = 2

% Ax+b<0 and x>=0 defines the feasible region
c = 10;
b1 = -4;
b2 = 5;
a11 = (c-b1)/c;
a21 = (c-b2)/c;

A = zeros(2,2);
b = zeros(2,1);

A(1,1) = a11;
A(1,2) = -1;
A(2,1) = -a21;
A(2,2) = 1;

b(1,1) = b1;
b(2,1) = -b2;

funArgs = [A,b];

% number of iterations
for j=1:1
    
    disp(['Dimenison: ',num2str(N), ', Iteration: ',num2str(j)])
    
    % Design centering mode
    inopts.mode = 0;
    
    inopts.r = 0.5;
    inopts.N_m = N*1e1;
    inopts.N_C = N*1e2;
    inopts.valP = 1/exp(1);
    inopts.c_T = 0.5;
    
    % Set-up of GaA options
    inopts.VerboseModulo = 2e3;
    inopts.bSaving='on';
    inopts.SavingModulo=1;
    inopts.bSaveCov = 0;
    
    inopts.LBounds = -2*ones(N,1);
    inopts.UBounds = 12*ones(N,1);
    inopts.MaxIter = N*1e4;
    
    inopts.Display = 'on';
    inopts.Plotting = 'on';
    
    inopts.funArgs = funArgs;
    
    % Start value, any feasible point in this case [1;1]
    xstart = [1;1];
    
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
end


% Latest covariance matrix
C = out.lastR^2*out.lastQ'*out.lastQ;

% Summary Plot
figure
[h,x,y,z]=error_ellipse('C',C(1:2,1:2),'mu',out.lastMu(1:2),'style','r');
set(h,'LineWidth',5)
hold on
plot(out.xAcc(:,1),out.xAcc(:,2),'b.','MarkerSize',0.4)
plot(out.lastMu(1,1),out.lastMu(2,1),'r.', 'MarkerSize',20)
grid on
xlabel('x_1')
ylabel('x_2')
title('Feasible points (blue) and GaA mean and covariance (red)')


% Compare true area with ellipsoidal area approximation
trueArea = b2*c + 1/2*(c-b2)*c - 1/2*(c+sign(b1)*c/(c-b1)*abs(b1))*c

% Scaling according to chi2 distribution
c_p=chi2inv(out.P_emp,N);

% Do an eigendecomposition of the true covariance
[eigVecs,eigVals] = eig(C);

% Rescale the covariance
elliAxes =  sqrt(c_p*eigVals);

% Compute the approximate volume of the minimum volume ellipsoid that
% corresponds to the statistical description
approxArea=prod(diag(elliAxes))*sphereVol(1,N)





