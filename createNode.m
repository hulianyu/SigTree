function node = createNode(objs_ID, min_pval, best_mcat, left, right)
    % Create a new tree node containing min_pval, left child, and right child
    node.objs = objs_ID;
    node.pval = min_pval;
    node.category = best_mcat;
    node.left = left;
    node.right = right;
end