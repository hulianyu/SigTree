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
load('CUBT_maxTree_Depth.mat')
load('CUBT_Ham_pi.mat')
load('CUBT_Ham_times.mat')
load('CUBT_Ham_Depth.mat')
CUBT_Ham_pi = struct2cell(CUBT_Ham_pi);
CUBT_Ham_metrics = zeros(18,2);
load('CUBT_MI_pi.mat')
load('CUBT_MI_times.mat')
load('CUBT_MI_Depth.mat')
CUBT_MI_pi = struct2cell(CUBT_MI_pi);
CUBT_MI_metrics = zeros(18,2);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    GT = X_data(:,1); %Ground Truth    
    %%
    results_Ham = ClusteringMeasure(GT, CUBT_Ham_pi{I,1});
    CUBT_Ham_metrics(I,1) = results_Ham(3); % Purity
    CUBT_Ham_metrics(I,2) = results_Ham(7); % F-score
    %%
    results_MI = ClusteringMeasure(GT, CUBT_MI_pi{I,1});
    CUBT_MI_metrics(I,1) = results_MI(3); % Purity
    CUBT_MI_metrics(I,2) = results_MI(7); % F-score
end