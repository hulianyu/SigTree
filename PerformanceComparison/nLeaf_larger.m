load('Depths_nLeaf.mat')
% For one-sided signrank
alpha = 0.05;
algo_nLeaf = {'K','SigTree','CUBT_{Max}','CUBT^{Ham}','CUBT^{MI}'};

t = Depths_nLeaf(:,2);
Depths_nLeaf(:,2) = Depths_nLeaf(:,1);
Depths_nLeaf(:,1) = t;
%% Larger
h_large = zeros(4,1);
for j = 2:5
    h_large(j-1,1) = signrank(Depths_nLeaf(:,j), Depths_nLeaf(:,1), 'Tail', 'right');
end
larger_algo = {'CUBT-Max'};

%% smaller
h_small = zeros(4,1);
for j = 2:5
    h_small(j-1,1) = signrank(Depths_nLeaf(:,j), Depths_nLeaf(:,1), 'Tail', 'left');
end
smaller_algo = [];