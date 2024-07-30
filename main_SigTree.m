addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/Evaluation']);
%% Load Data sets
% filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
%     'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
%     'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
%     'mammographic-mass','tic-tac-toe','car','chess','mushroom','nursery','connect4','splice','autism',...
%     'onlineshoppersintention','NPHAdoctorvisits');
% rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
%     'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce', 'Ch', 'Mu', 'Nu', 'Cf', 'Sp', 'Au',...
%     'Os', 'Np'};
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
IS = size(filename,1);
Results = zeros(IS,3);
eTime = zeros(IS,1);
numPi = zeros(IS,2);
Depths = zeros(IS,2);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    GT = X_data(:,1); %Ground Truth
    K = length(unique(GT));
    [N,M] = size(X);
    Q = 0;
    for m =1:M
        Q = Q + length(unique(X(:,m)));
    end
    objsID = 1:N;
    X = [objsID' X];
    num_node = 0;
    pi_Node = zeros(N+1,1);
    tic
    [Node, pi_Node] = Binary_divide(X,Q,pi_Node); 
    % [Node, pi_Node] = Binary_divide_K(X,Q,pi_Node,K); 
    eTime(I,1) = toc;
    pi = pi_Node(1:end-1);
    %% Clustering quality
    results = ClusteringMeasure(GT, pi);
    Results(I,1) = results(3); % Purity
    Results(I,2) = results(7); % F-score
    Results(I,3) = results(4); % ARI
    %% Depths
    numPi(I,1) = length(unique(GT));
    numPi(I,2) = length(unique(pi));
    Depths(I,1) = treeDepth(Node);
    Depths(I,2) = averageLeafDepth(Node);
end
%     drawTree(Node, rowNames, I); % Call the function to draw the tree