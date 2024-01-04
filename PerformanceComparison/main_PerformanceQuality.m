load('PerformanceQuality.mat')
% For Bonferroni-Dunn (BD) test
alpha = 0.05;

%% All algorithms (except for DV)
algo = {'SigTree','CUBT^{Ham}','CUBT^{MI}','k-modes','Entropy','CDE','CDCDR','IMM','RDM','SHA'};
Purity_list = [Metric_SigTree(:,1) CUBT_Ham_metrics(:,1) CUBT_MI_metrics(:,1) Metric_kmodes(:,3) Metric_Entropy(:,3) Metric_CDE(:,3) Metric_CDC_DR(:,3) IMM_metrics(:,1) RDM_metrics(:,1) SHA_metrics(:,1)];
Fscore_list = [Metric_SigTree(:,2) CUBT_Ham_metrics(:,2) CUBT_MI_metrics(:,2) Metric_kmodes(:,7) Metric_Entropy(:,7) Metric_CDE(:,7) Metric_CDC_DR(:,7) IMM_metrics(:,2) RDM_metrics(:,2) SHA_metrics(:,2)];
[cd,f] = criticaldifference('Purity_CD',Purity_list,algo,alpha);
[cd1,f1] = criticaldifference('Fscore_CD',Fscore_list,algo,alpha);

%% Interpretable clustering algorithms
algo = {'SigTree','CUBT^{Ham}','CUBT^{MI}','IMM','RDM','SHA'};
Purity_list_Tree = [Metric_SigTree(:,1) CUBT_Ham_metrics(:,1) CUBT_MI_metrics(:,1) IMM_metrics(:,1) RDM_metrics(:,1) SHA_metrics(:,1)];
Fscore_list_Tree = [Metric_SigTree(:,2) CUBT_Ham_metrics(:,2) CUBT_MI_metrics(:,2) IMM_metrics(:,2) RDM_metrics(:,2) SHA_metrics(:,2)];
% [cd,f] = criticaldifference('PurityTree_CD',Purity_list,algo,alpha);
[cd2,f2] = criticaldifference('FscoreTree_CD',Fscore_list_Tree,algo,alpha);