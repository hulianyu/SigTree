function [oc,chi] = chi_squared_test_table(X,pi)
[N,M] = size(X);
k_list = unique(pi);
K = length(k_list);
p = zeros(1,M);
chi = 0;
for m=1:M
    Xm = X(:,m);
    Q = length(unique(Xm));
    oc = zeros(Q+1,K+1);% observed counts
    for k=1:K
        for q=1:Q
            set = Xm(pi==k_list(k));
            oc(q,k) = sum(set==q);
        end
    end
    oc(Q+1,:) = sum(oc);
    oc(:,K+1) = sum(oc,2);
    ec = (oc(end,1:end-1).*oc(1:end-1,end))./N;
    oc = oc(1:end-1,1:end-1);
    s = sum(((oc-ec).^2)./ec,'all');% chi-squared statistic
    chi = chi + s;
    df = (Q-1)*(K-1);% degrees of freedom
%     p(1,m) = chi2cdf(s,df,'upper');
    p(1,m) = simplifiedChi2cdf(s,df); % 'upper'
end
end