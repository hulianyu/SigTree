function pi = CDC_DR_SC(data,k)
% [D,~]=CDC_DR(data,k,'AE','Joint');
[D,~]=CDC_DR(data,k,'SC','Joint');
pi=kmeans(D,k);
end

function [D,TS]=CDC_DR(data,k,Type,Type2)
data = trans2unique(data);
tic1 = tic;
[r,c]=size(data);
[X]=One_Hot(data);
q=fix(length(X(1,:))/c);
%q=k;
S=Similarity(X);
%q=k;
switch Type
    case 'CDE'
        D=CDE(data);
    case 'NMF'
      [H,W]=nnmf(S,k,'replicates',50);
      %H=W;
      H= bsxfun(@rdivide,H,sqrt(sum(H.^2,2)) + 1e-10);
      D=Pooling(X,H,Type2);
    case 'KM'
      [~,H] = kmeans(S,k,'Distance','cosine','EmptyAction','singleton','Replicates',10); 
      H=H';
      H= bsxfun(@rdivide,H,sqrt(sum(H.^2,2)) + 1e-10);
      D=Pooling(X,H,Type2);
    case 'SC'
       %S=S-diag(diag(S));
       degs = sum(S, 2); 
       D=diag(degs);
       %L=D-S;
       degs(degs == 0) = 1e-10;
       D = full(spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2)));
       S = D * S * D;
       %diff=1e-10;
%        [H, ~] = eigs(L, q, diff);
       
       [A, B] = eig(S);
       [~,R]=sort(diag(B),'descend');
       H=A(:,R(1:k));
       %H=A;
       H= bsxfun(@rdivide,H,sqrt(sum(H.^2,2)) + 1e-10);
       D=Pooling(X,H,Type2); 
    case 'DI'
       H=S;
       H= bsxfun(@rdivide,H,sqrt(sum(H.^2,2)) + 1e-10);
       D=Pooling(X,H,Type2);
    case 'AE'
       autoenc =trainAutoencoder(S,q,'MaxEpochs',1000);
       H = encode(autoenc,S);
       
       %autoenc1 =trainAutoencoder(H,fix(q/2),'MaxEpochs',1000);
       %stacknet = stack(autoenc,autoenc1);
       %deepnet = train(stacknet,S);
       %H=stacknet(S);
       H=H';
       H= bsxfun(@rdivide,H,sqrt(sum(H.^2,2)) + 1e-10);
       D=Pooling(X,H,Type2);
end
%I = kmeans(D,k,'Replicates',10);
%[~,I]=max(D,[],2);
TS=toc(tic1);
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

function S=Similarity(X)
Q=X';
S=1-squareform(pdist(Q,'cosine' ));
%S=squareform(exp(-pdist(Q,'squaredeuclidean')/0.01));
% [V,R]=sort(S,2,'descend');
% [r,~]=size(S);
% Sk=zeros(r,r);
% for i=1:r
%     Sk(i,R(i,[1:k]))=S(i,R(i,[1:k]));
% end
end

function D=Pooling(X,H,Type)
switch Type
    case 'Mean'
        D=X*H;
    case 'Joint'
        [n,~]=size(X);
        for i=1:n
            I=find(X(i,:)>0);
            T=[];
            for j=1:length(I)
                T=[T,H(I(j),:)];
            end
            D(i,:)=T;
        end
end
end

function [X]=One_Hot(E)
[r,cc]=size(E);
Kh=zeros(1,cc);
Ks=zeros(1,cc+1);
for i=1:cc
    cl{i}=unique(E(:,i));
    Kh(i)=length(cl{i});
    Ks(i+1)=Ks(i)+Kh(i);
end
nc=sum(Kh);
X=zeros(r,nc);
for j=1:cc
    for i=1:r
        X(i,Ks(j)+find(cl{j}==E(i,j)))=1;
    end
end
end
