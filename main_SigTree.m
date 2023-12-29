addpath([cd '/']);
addpath([cd '/Datasets']);
addpath([cd '/Evaluation']);
%% Load Data sets
filename = char('lenses','lung-cancer','soybean-small','zoo','dna-promoter',...
    'hayes-roth','lymphography','heart-disease','solar-flare','primary-tumor',...
    'dermatology','house-votes','balance-scale','credit-approval','breast-cancer-wisconsin',...
    'mammographic-mass','tic-tac-toe','car');
rowNames = {'Ls', 'Lc', 'So', 'Zo', 'Ps', 'Hr', 'Ly', 'Hd', 'Sf',...
    'Pt', 'De', 'Hv', 'Bs', 'Ca', 'Bc', 'Mm', 'Tt', 'Ce'};
Results = zeros(18,2);
eTime = zeros(18,1);
numPi = zeros(18,2);
Depths = zeros(18,2);
for I=1:18
    disp(I);
    X_data = load([strtrim(filename(I,:)), '.txt']); %Load a Dataset
    X = X_data(:,2:end); %Data set
    GT = X_data(:,1); %Ground Truth 
    [N,M] = size(X);
    Q = 0;
    for m =1:M
        Q = Q + length(unique(X(:,m)));
    end
    objsID = 1:N;
    X = [objsID' X];
    tic
    Node = Binary_divide(X,Q);
    if isempty(Node.category)
        T = 1;
        Node = Binary_divide_nonsig(X,T);
    end
    pi_Node = AssignLeaf(Node,N);
    eTime(I,1) = toc;
    %% Clustering quality
    results = ClusteringMeasure(GT, pi_Node);
    Results(I,1) = results(3); % Purity
    Results(I,2) = results(7); % F-score
    %% Depths
    numPi(I,1) = length(unique(GT));
    numPi(I,2) = length(unique(pi_Node));
    Depths(I,1) = treeDepth(Node);
    Depths(I,2) = averageLeafDepth(Node);
end
%     drawTree(Node, rowNames, I); % Call the function to draw the tree