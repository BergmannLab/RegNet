% Script that runs Gaussian Adaptation on various target distributions
clear inopts
close all

% Example from Liang et al.: a 2D 20-funnel distribution
% funName = 'fLiangEx1'
% global LiangExMat
% load LiangExMat.mat

% Example from Liang et al.: a nD double-funnel distribution
% funName = 'fLiangEx2'

% Example from Haario et al., 1999: a 8D banana-shaped distribution
funName = 'fHaario';

% Haario's banana-shaped distribution needs: 
% funArgs: 1 (untwisted), 3(moderatly twisted) or 4(strongly twisted)
funArgs = 3;

i=1;

% Dimension
N = 8

numIter = 10;
MaxSample = 4e4;
sampleData = zeros(MaxSample,N,numIter); 

tempSum=0;
% number of iterations
for j=1:numIter
    j
    
    % MCMC sampler
    inopts.mode = 2;
    
    % Set-up of GaA options
    inopts.VerboseModulo = Inf;
    inopts.bSaving='on';
    inopts.SavingModulo=1;
    inopts.LBounds = -80*ones(N,1);
    inopts.UBounds = 80*ones(N,1);
    inopts.r = 1;
    inopts.N_mu = 1;
    inopts.MaxIter = MaxSample;
    
    inopts.valP = 0.234;
    inopts.TolFun = 0;
    inopts.TolX = 0;
    
    % Comment out if Liang's target functions are used
    inopts.funArgs = funArgs;
    
    inopts.Display = 'off';
    inopts.Plotting = 'off';
    
    % Start value (as in the CEC 2010 paper)
    xstart = -1 + 2*rand(N,1);
    
    [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);
    sampleData(:,:,j)=out.xRaw;
    
end

% Plot a 2D projection of a randomly selected sample of the distribution
figure

randInd = ceil(rand*numIter)
currData = squeeze(sampleData(:,:,randInd));

plot(currData(:,1),currData(:,2),'b.','MarkerSize',0.5)
hold on
grid on
xlabel('x_1')
ylabel('x_2')
title('Randomly selected sample from the target distribution')






