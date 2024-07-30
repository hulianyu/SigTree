function [node, pi] = Binary_divide(X,Q,pi)
node = createNode(X(:,1), [], [], [], []); % Create a new node
pi(end) = pi(end)+1; % Record the cluster number
pi(node.objs) = pi(end);  % Assign a leaf to a cluster
% disp(pi(end));
if  size(X,1)<=5
    return;
end
% Calculate min_pval and best partition index
[min_pval,best_pi,best_m,best_cat] = find_best_mq(X);
% If min_pval is greater than the threshold
if pi(end) ~= 1 && min_pval>0.01/Q
    return;
end
%%
% Recursively apply Binary_divide to the left and right sub-arrays
if min(simplifiedHistcounts(best_pi, 0.5:2.5))>=3
    X_left = X(best_pi == 1,:);
    X_right = X(best_pi == 2,:);
    node.pval = min_pval;
    node.category = [best_m best_cat];
    [node.left, pi] = Binary_divide(X_left,Q,pi);
    [node.right, pi] = Binary_divide(X_right,Q,pi);
end
end