load('RunningTimes_list.mat')
% For one-sided signrank
% alpha = 0.05;
algo = {'SigTree','CUBT^{Ham}','CUBT^{MI}','k-modes','Entropy','CDE','CDC\_DR','IMM','RDM','SHA'};

%% Faster
h_fast = zeros(10,1);
for j = 2:11
    h_fast(j-1,1) = signrank(RunningTimes_list(:,j), RunningTimes_list(:,1), 'Tail', 'left');
end
Faster_algo = {'CDE','CDC\_DR','IMM','RDM','SHA'};
%% Slower
h_slow = zeros(10,1);
for j = 2:11
    h_slow(j-1,1) = signrank(RunningTimes_list(:, 1), RunningTimes_list(:, j),'tail', 'left');
end
Slower_algo = {'DV','CUBT^{Ham}','CUBT^{MI}','Entropy'};
