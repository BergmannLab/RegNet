addpath('./matlab/optim/GaA/');

ub = 10; % upper bound
lb = -10;% lower bound
sample = {};
for i=1:npar,
RandStream.setGlobalStream(RandStream.create('mlfg6331_64','seed',rand+1000000*sum(clock)+i));
opts={};
opts.MaxIter = 70000;
opts.VerboseModulo = 10000;
opts.r = 1;
opts.mode = 1;

sockname =  [socketname '_' num2str(i)];
opts.funArgs = sockname;  

dim = getDimension(sockname);

opts.LBounds       = lb *ones(1,dim)';   
opts.UBounds       = ub *ones(1,dim)'; 
opts.Plotting = 'off';
opts.bRestart = 0;
par  =-2 *2*rand(1,dim);
[xmin,fmin,coutval,out]= gaussAdapt('requestEval',par',opts);
if(fmin<fitok)
    fnamei=[fname num2str(i)];
   opts.StopFitness = -Inf;
   opts.mode = 2;
   fnamei=[fname num2str(i) '.par'];
   fd = fopen(fnamei,'w');
   if(fast==0)
          opts.MaxIter = 10000;
   opts.SavingModulo  = 1;
   [xmin,fmin,coutval,out]= gaussAdapt('requestEvalMCMC',xmin,opts);
   sample{i} = out.xAcc;
   fwrite(fd,out.xAcc','float32');
   fclose(fd);
    dlmwrite([fname num2str(i) '.map'],[fmin; xmin]','precision',10);
   else
       sample{i} = xmin;
       fwrite(fd,xmin','float32');
       fclose(fd);
       closeServer(sockname);
       break;
   end
   end
   closeServer(sockname);
end
files =[];
mapfiles=[];
if(length(sample)>0),
for i=1:size(sample,2),
    if(length(sample{i})>0)
        files = [files fname num2str(i) '.par '];
        mapfiles = [mapfiles fname num2str(i) '.map '];
    end
end
unix(['cat ' files '> ' fname '.par']);
unix(['rm -f ' files]);
unix(['cat ' mapfiles '> ' fname '_map.txt']);
unix(['rm -f ' mapfiles]);
end
