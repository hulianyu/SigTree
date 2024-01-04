load('Depths_list.mat')
% For one-sided signrank
alpha = 0.05;
algo_avgDepth = {'SigTree','CUBT_{Max}','CUBT^{Ham}','CUBT^{MI}','IMM','RDM','SHA'};

% #[1]SHA [2]RDM [3]IMM
avgDepth_list = [SigTree_Depth(:,3) CUBT_maxTree_Depth(:,3) CUBT_Ham_Depth(:,3) CUBT_MI_Depth(:,3) AvgLeafDepth(:,3) AvgLeafDepth(:,2) AvgLeafDepth(:,1)];
%% Larger
avgDepth_large = zeros(6,1);
for j = 2:7
    avgDepth_large(j-1,1) = signrank(avgDepth_list(:,j), avgDepth_list(:,1), 'Tail', 'right');
end
larger_algo = {'CUBT-Max','CUBT^{MI}'};

%% smaller
avgDepth_small = zeros(6,1);
for j = 2:7
    avgDepth_small(j-1,1) = signrank(avgDepth_list(:,j), avgDepth_list(:,1), 'Tail', 'left');
end
smaller_algo = [];