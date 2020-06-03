% General settings
global initial_flag

clear inopts
close all

% Function
funName = 'benchmark_func'

% Boundaries
cecPath = 'CECPATH/cec2005/';
load([cecPath,'benchmarkBounds.mat'])

% Minima
minimaVec=load([cecPath,'fbias_data']);
% Accuracy
accVec=load([cecPath,'acc_data']);

% max function number
maxFunc=15;

stopFitness = minimaVec.f_bias(1:maxFunc)' + accVec.acc_data(1:maxFunc);                                                                                                                                                                         
                                                                                                                                                                                                                                                 
% Save data                                                                                                                                                                                                                                      
outCEC2005=cell(maxFunc,25);                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                 
% Dimension                                                                                                                                                                                                                                      
N = 10                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                 
    tic                                                                                                                                                                                                                                          
    % Function to be optimized                                                                                                                                                                                                                   
for i=1:maxFunc                                                                                                                                                                                                            
                                                                                                                                                                                                                                                 
        initial_flag=0                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                 
        for j=1:25                                                                                                                                                                                                                               
          disp(['Current state: ',num2str(i),', ',num2str(j)])                                                                                                                                                                                   
          % Set-up of GaA options                                                                                                                                                                                                                
          inopts.VerboseModulo = 2*N*100;                                                                                                                                                                                                            
          inopts.bSaving='on';                                                                                                                                                                                                                   
          inopts.SavingModulo=N*100;                                                                                                                                                                                                              
          inopts.mode = 1;
          inopts.LBounds = benchmarkBounds(i,1)*ones(N,1);                                                                                                                                                                                                 
          inopts.UBounds = benchmarkBounds(i,2)*ones(N,1);
          inopts.BoundPenalty=0;
          inopts.r = 1/exp(1)*max((inopts.UBounds-inopts.LBounds));
          inopts.MaxIter = 1e4*N;                                                                                                                                                                                                                
          inopts.StopFitness = stopFitness(i);                                                                                                                                                                                                   
          inopts.valP = 1/exp(1);                                                                                                                                                                                                                
          inopts.Display = 'off';                                                                                                                                                                                                                
          inopts.Plotting = 'off';  
          inopts.funArgs = i;
          inopts.ThreshRank = 0;
          inopts.N_T = exp(1)*N;
          inopts.bRestart = 1;
          inopts.inc_T = 1;
          
                    
          % Start value                                                                                                                                                                                                                          
          xstart = inopts.LBounds + rand(N,1).*(inopts.UBounds-inopts.LBounds);
          % Special function                                                                                                                                                                                                                     
          if(i==7)               
              inopts.UBounds=600*ones(N,1);                                                                                                                                                                                                              
              inopts.LBounds=-600*ones(N,1);                                                                                                                                                                                                               
          end                                                                                                                                                                                                                                    
          [xmin,fmin,counteval,out] = gaussAdapt(funName,xstart,inopts);                                                                                                                                                                         
          outCEC2005{i,j} = out;                                                                                                                                                                                                                 
        end                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                 
    end                                                                                                                                                                                                                                          
    toc                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                 
%save('outCEC2005','outCEC2005')         