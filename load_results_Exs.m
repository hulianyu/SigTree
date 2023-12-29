addpath([cd '/']);
addpath([cd '/ComparedMethods']);
addpath([cd '/Datasets']);
addpath([cd '/Evaluation']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
cfilename = char('lenses','lung_cancer','soybean_small','zoo','dna_promoter',...
    'hayes_roth','lymphography','heart_disease','solar_flare','primary_tumor',...
    'dermatology','house_votes','balance_scale','credit_approval','breast_cancer_wisconsin',...
    'mammographic_mass','tic_tac_toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
%% Load Results
load('SHA_pi.mat')
load('RDM_pi.mat')
load('IMM_pi.mat')
SHA_metrics = zeros(18,2);
RDM_metrics = zeros(18,2);
IMM_metrics = zeros(18,2);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    GT = X_data(:,1); %Ground Truth
    for run=1:50
        j = run + (I-1)*50;
        %%
        SHA_pi{1,j} = (SHA_pi{1,j}+1)';
        results_SHA = ClusteringMeasure(GT, SHA_pi{1,j});
        SHA_metrics(I,1) = SHA_metrics(I,1) + results_SHA(3); % Purity
        SHA_metrics(I,2) = SHA_metrics(I,2) + results_SHA(7); % F-score
        %%
        RDM_pi{1,j} = (RDM_pi{1,j}+1)';
        results_RDM = ClusteringMeasure(GT, RDM_pi{1,j});
        RDM_metrics(I,1) = RDM_metrics(I,1) + results_RDM(3); % Purity
        RDM_metrics(I,2) = RDM_metrics(I,2) + results_RDM(7); % F-score
        %%
        IMM_pi{1,j} = (IMM_pi{1,j}+1)';
        results_IMM = ClusteringMeasure(GT, IMM_pi{1,j});
        IMM_metrics(I,1) = IMM_metrics(I,1) + results_IMM(3); % Purity
        IMM_metrics(I,2) = IMM_metrics(I,2) + results_IMM(7); % F-score
    end
end
%%
SHA_metrics = SHA_metrics/50;
RDM_metrics = RDM_metrics/50;
IMM_metrics = IMM_metrics/50;
disp([mean(SHA_metrics,1);mean(RDM_metrics,1);mean(IMM_metrics,1)])
load('Execution_times.mat')
load('NumLeaf.mat')
load('MaxDepth.mat')
load('AvgLeafDepth.mat')