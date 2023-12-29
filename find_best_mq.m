function [minpval,bpi,bm,bcat] = find_best_mq(X)
X(:,1) = []; % objsID
M = size(X,2);
pv = cell(M,1);
pi = cell(M,1);
discat = cell(M,1);
for m=1:M
    % keep the order in discat
    [discat{m,1}, ~, X(:,m)] = unique(X(:,m), 'stable');
end
for m=1:M
    Q = length(discat{m,1});
    if Q~=1
        for q=1:Q
            [pi_current, pv_current] = sum_chi(X,m,q); % add Yates correction
            pi{m,1} = [pi{m,1} pi_current];
            pv{m,1} = [pv{m,1} pv_current];
        end
    else
        pi{m,1} = ones(size(X,1),1);
        pv{m,1} = 1;
    end
end
% find the best m
[minpval, bm] = min(cellfun(@min, pv));
% find the best category according to the recorded discat
bq = find(pv{bm}==minpval,1);
bcat = discat{bm}(bq,1);
% find the best partition
bpi = pi{bm}(:,bq);
end