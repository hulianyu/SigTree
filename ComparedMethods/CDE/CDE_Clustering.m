function pi = CDE_Clustering(X,K)
X_u = trans2unique(X);
new_rep=CDE(X_u);
pi=kmeans(new_rep,K);
end

function X_u = trans2unique(X)
M = size(X,2);
X_u = X;
Q = max(X(:,1:M));
for m=2:M
    Qm = sum(Q(1:m-1));
    X_u(:,m) = X(:,m) + Qm;
end
end