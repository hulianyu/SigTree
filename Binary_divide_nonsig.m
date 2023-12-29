function node = Binary_divide_nonsig(X,T)
node = createNode(X(:,1), [], [], [], []); % Create a new node
if  size(X,1)<=5
    return;
end
% Calculate min_pval and best partition index
[min_pval,best_pi,best_m,best_cat] = find_best_mq(X);
% If min_pval is greater than the threshold
if  T > 1
    return;
end
T = 2;
%%
% Recursively apply Binary_divide to the left and right sub-arrays
X_left = X(best_pi == 1,:);
X_right = X(best_pi == 2,:);
if height(X_left)>2 && height(X_right)>2
    node.pval = min_pval;
    node.category = [best_m best_cat];
    node.left = Binary_divide_nonsig(X_left,T);
    node.right = Binary_divide_nonsig(X_right,T);
end
end