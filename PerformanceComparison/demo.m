load('PerformanceQuality.mat')
% For Bonferroni-Dunn (BD) test
alpha = 0.1;
algo = {'SigTree','CUBT^{Ham}','CUBT^{MI}','k-modes','Entropy','CDE','CDC\_DR','IMM','RDM','SHA'};
Purity_list = [Metric_SigTree(:,1) CUBT_Ham_metrics(:,1) CUBT_MI_metrics(:,1) Metric_kmodes(:,3) Metric_Entropy(:,3) Metric_CDE(:,3) Metric_CDC_DR(:,3) IMM_metrics(:,1) RDM_metrics(:,1) SHA_metrics(:,1)];
Fscore_list = [Metric_SigTree(:,2) CUBT_Ham_metrics(:,2) CUBT_MI_metrics(:,2) Metric_kmodes(:,7) Metric_Entropy(:,7) Metric_CDE(:,7) Metric_CDC_DR(:,7) IMM_metrics(:,2) RDM_metrics(:,2) SHA_metrics(:,2)];
[cd,f] = criticaldifference('Purity_CD',Purity_list,algo,alpha);
[cd,f] = criticaldifference('Fscore_CD',Fscore_list,algo,alpha);

% %% Significant (DV) data sets
% list = [3 4 7 8 11 12 15];
% algo = {'SigTree','DV','CUBT^{Ham}','CUBT^{MI}'};
% Purity_list = [Metric_SigTree(list,1) Metric_DV(list,3) CUBT_Ham_metrics(list,1) CUBT_MI_metrics(list,1)];
% Fscore_list = [Metric_SigTree(list,2) Metric_DV(list,7) CUBT_Ham_metrics(list,2) CUBT_MI_metrics(list,2)];
% [cd,f] = criticaldifference('PuritySig_CD',Purity_list,algo,alpha);
% [cd,f] = criticaldifference('FscoreSig_CD',Fscore_list,algo,alpha);

algo = {'SigTree','CUBT^{Ham}','CUBT^{MI}','IMM','RDM','SHA'};
Purity_list_Tree = [Metric_SigTree(:,1) CUBT_Ham_metrics(:,1) CUBT_MI_metrics(:,1) IMM_metrics(:,1) RDM_metrics(:,1) SHA_metrics(:,1)];
Fscore_list_Tree = [Metric_SigTree(:,2) CUBT_Ham_metrics(:,2) CUBT_MI_metrics(:,2) IMM_metrics(:,2) RDM_metrics(:,2) SHA_metrics(:,2)];
% [cd,f] = criticaldifference('PurityTree_CD',Purity_list,algo,alpha);
[cd,f] = criticaldifference('FscoreTree_CD',Fscore_list_Tree,algo,alpha);