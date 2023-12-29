% Calculate the depth of the tree
function maxDepth = treeDepth(node)
if isempty(node) || (~isstruct(node.left) && ~isstruct(node.right))
    maxDepth = 0;
else
    maxDepthLeft = treeDepth(node.left);
    maxDepthRight = treeDepth(node.right);
    maxDepth = 1 + max(maxDepthLeft, maxDepthRight);
end
end
