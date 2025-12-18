%%%%%%%%%%%%%%%%%%%%%%
% Export HMM results %
%%%%%%%%%%%%%%%%%%%%%%

%set parameters
K =12; % number of states
T = Ti; % Max time length
rep = 1; % repetition of minimum free energy

load(['..\.\results\HMMrun_rep_' num2str(rep) '.mat'])

cd("..\.\features")

%extract features
ave_actives = getFractionalOccupancy(Gamma,T,options);

getStateIntervalTimes(vpath,T,options);

swithchs = getSwitchingRate(Gamma,T,hmm.train);

each_pop{length(T),K} = [];

each_LT = getStateLifeTimes(vpath,Ti,options,[],[],0 );

for i = 1:length(T)
    for j = 1:K
        each_pop{i,j} = size(each_LT{i,j},2);
    end
end

mean_each_LT = cellfun(@mean,each_LT);

%transition matrix
accum{length(T),1} = [];
accum{1,1} = sum(T{1});
for i = 2:length(T)
    accum{i,1} = accum{i-1,1} + sum(T{i});
end

M{length(T),1} = [];
M{1,1} = 1:T{1};
for i = 2:length(T)
    M{i,1} = (accum{i-1}+1):accum{i};
end

[Gamma, Xi] = hmmdecode(f1, T, hmm, 0);
Trans_probs = getMaskedTransProbMats(cell2mat(f1), cell2mat(T'),hmm,M,Gamma,Xi);

mkdir tp_raw
cd(".\tp_raw")
t_data = Trans_probs';
for i = 1:length(t_data)
    writematrix(t_data{i,1}, strcat("tmat", num2str(i,"%02d"), "_.csv"));
end

%save features
folderInfo = dir;
folderInfo = folderInfo(~ismember({folderInfo.name},{'.','..'}));
folderlist = {folderInfo(1:length(folderInfo)).name};

N=length(folderlist);
demo_cell = cell(N,1);

for i = 1:N
  mats =  readmatrix(folderlist{i}); 
  demo_cell{i,1} = reshape(mats,[1,K^2]);
end

resmat = cell2mat(demo_cell);

cd("..\.\features")

writematrix(ave_actives, "ave_actives.csv");
writematrix(mean_each_LT, "each_LT.csv");
writecell(each_pop, "each_pop.csv");
writematrix(swithchs, "switches.csv");
writematrix(resmat, "transprobs.csv");







