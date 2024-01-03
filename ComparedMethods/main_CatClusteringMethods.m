clearvars
addpath([cd '/']);
%% 18 UCI data sets
addpath([cd '/Datasets']);
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
IS = size(filename,1);
RunningTimes = zeros(IS,5); % seconds (for each clustering algorithm)
Pi_random = cell(IS,1);
%% Load Evaluation package
addpath([cd '/Evaluation']);
%% set execute times (for each clustering algorithm)
ET = 50;
%% Comparison methods: k-modes, Entropy, CDE, CDC_DR, DV
% [1] k-modes
addpath([cd '/ComparedMethods/k-modes']);
Pi_kmodes = cell(IS,1);
Metric_kmodes = zeros(IS,7);
% [2] Entropy
addpath([cd '/ComparedMethods/Entropy']);
Pi_Entropy = cell(IS,1);
Metric_Entropy = zeros(IS,7);
% [3] CDE
addpath([cd '/ComparedMethods/CDE']);
Pi_CDE = cell(IS,1);
Metric_CDE = zeros(IS,7);
% [4] CDC_DR
addpath([cd '/ComparedMethods/CDC_DR']);
Pi_CDC_DR = cell(IS,1);
Metric_CDC_DR = zeros(IS,7);
% [5] DV
addpath([cd '/ComparedMethods/DV']);
Pi_DV = cell(IS,1);
Metric_DV = zeros(IS,7);
DV_k = zeros(IS,1); % Detect k clusters
%% Choose a data set I
for I=4:4
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    [N,M] = size(X);
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT)); %Cluster Number
    %% Performance
    %% [1] k-modes
    disp("k-modes");
    Pi = zeros(size(X,1),ET);
    tic
    for runs = 1:ET
        pi_runs = kmode(X, K);
        Pi(:,runs) = pi_runs;
    end
    RunningTimes(I,1) = toc;
    RunningTimes(I,1) = RunningTimes(I,1)/ET;
    Pi_kmodes{I,1} = Pi;
    % evaluate the metrics (average)
    for runs = 1:ET
        Metric_kmodes(I,:) = Metric_kmodes(I,:) + [ClusteringMeasure(GT, Pi(:,runs))];
    end
    Metric_kmodes(I,:) = Metric_kmodes(I,:)./ET;
    %% [2] Entropy
    disp("Entropy");
    % initialization of pi
    Pi = zeros(size(X,1),ET);
    tic
    [eva0,freq_pi0,freq_m0,eva_m0,eva_2_0] = EE_init(X,K);
    for runs = 1:ET
        [~,pi_runs] = minMC_EE(X,K,eva0,freq_pi0,freq_m0,eva_m0,eva_2_0);
        Pi(:,runs) = pi_runs;
    end
    RunningTimes(I,2) = toc;
    RunningTimes(I,2) = RunningTimes(I,2)/ET;
    Pi_Entropy{I,1} = Pi;
    % evaluate the metrics (average)
    for runs = 1:ET
        Metric_Entropy(I,:) = Metric_Entropy(I,:) + [ClusteringMeasure(GT, Pi(:,runs))];
    end
    Metric_Entropy(I,:) = Metric_Entropy(I,:)./ET;
    %% [3] CDE
    disp("CDE");
    Pi = zeros(size(X,1),ET);
    tic
    for runs = 1:ET
        pi_runs = CDE_Clustering(X,K);
        Pi(:,runs) = pi_runs;
    end
    RunningTimes(I,3) = toc;
    RunningTimes(I,3) = RunningTimes(I,3)/ET;
    Pi_CDE{I,1} = Pi;
    % evaluate the metrics (average)
    for runs = 1:ET
        Metric_CDE(I,:) = Metric_CDE(I,:) + [ClusteringMeasure(GT, Pi(:,runs))];
    end
    Metric_CDE(I,:) = Metric_CDE(I,:)./ET;
    %% [4] CDC_DR
    disp("CDC_DR");
    Pi = zeros(size(X,1),ET);
    tic
    for runs = 1:ET
        pi_runs = CDC_DR_SC(X,K);
        Pi(:,runs) = pi_runs;
    end
    RunningTimes(I,4) = toc;
    RunningTimes(I,4) = RunningTimes(I,4)/ET;
    Pi_CDC_DR{I,1} = Pi;
    % evaluate the metrics (average)
    for runs = 1:ET
        Metric_CDC_DR(I,:) = Metric_CDC_DR(I,:) + [ClusteringMeasure(GT, Pi(:,runs))];
    end
    Metric_CDC_DR(I,:) = Metric_CDC_DR(I,:)./ET;
    %% [5] DV
    disp("DV");
    tic
    try
        Pi = ccdv(X);
    catch
        Pi = 0; % fail to detect clusters
    end
    RunningTimes(I,5) = toc;
    if Pi~=0
        Pi_DV{I,1} = Pi;
        % evaluate the metrics
        Metric_DV(I,:) = ClusteringMeasure(GT, Pi);
    else
        Metric_DV(I,:) = 0;
    end
    DV_k(I,1) = max(Pi);
    %% Save the results
%     save('_CatClusteringResults.mat','RunningTimes','ET','Pi_kmodes','Metric_kmodes',...
%         'Pi_Entropy','Metric_Entropy','Pi_CDC_DR','Metric_CDC_DR','Pi_CDE','Metric_CDE',...
%         'Pi_DV','Metric_DV','DV_k');
end