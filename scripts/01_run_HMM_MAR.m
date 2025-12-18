%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse resting-state fMRI data using HMM-MAR %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup specific parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%directory and others
DirIn = '..\.\data'; %data directory
DirOut = '..\.\results'; %output directory

repetitions = 150; % to run estimation multiple times 
TR = 0.8;  

% options
options = struct();
options.order = 0; % no autoregressive components
options.zeromean = 0; % model the mean
options.covtype = 'uniquefull'; % full covariance matrix
options.Fs = 1/TR; %sampling frequency
options.verbose = 1; %show program process
options.standardise = 1;
options.inittype = 'HMM-MAR';
options.tol = 1e-5; 
options.cyc = 1000; %default = 1000
options.initcyc = 100; %default = 100
options.initrep = 5; %default = 5
options.leida = 0; %leading eingenvector dynamics analysis 
options.useParallel = 0; 
options.pca = .90;

% stochastic options
options.BIGNbatch = 100;
options.BIGtol = 1e-5;
options.BIGcyc = 1000;
options.BIGundertol_tostop = 10;
options.BIGbase_weights = 0.9;

%%%%%%%%%%%%%%%
% data import %
%%%%%%%%%%%%%%%

cd(DirIn)
folderdata = dir;
folderdata = folderdata(~ismember({folderdata.name},{'.','..'}));
datafolderlist = {folderdata.name};
[c, subN] = size(datafolderlist);

f1{subN,1} = {};
T1{subN,1} = {};

for i = 1:subN
    f1{i,1} = readmatrix(datafolderlist{1,i}, 'Range',[2 2]);
    Ti{i,1} =  size(f1{i,1},1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct estimation multiple times for minimum free envergy %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fe0 = cell(20,1);

for s = 1:20
    try
    options.K = s; % number of states 
    [hmm, Gamma, ~, vpath] = hmmmar(f1,Ti,options);
    fe0{s,1} = hmmfe(f1,Ti,hmm,Gamma);
    catch 
    end
end

[data_min0, idx_min0] = min(cell2mat(fe0));
idx_min0

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate HMM parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.K = idx_min; % set number of states 

for r = 1:repetitions
    try
    [hmm, Gamma, ~, vpath] = hmmmar(f1,Ti,options);
    save([DirOut 'HMMrun_rep_' num2str(r) '.mat'],'Gamma','vpath','hmm');
    fe1{r,1} = hmmfe(f1,Ti,hmm,Gamma);
    catch 
    end
end

%select minimum free energy model
[data_min1, idx_min1] = min(cell2mat(fe1));
idx_min1


