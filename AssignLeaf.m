function pi = AssignLeaf(node,N)
pi = zeros(N,1);
Clusters = {};  % Initialize the cell array
Clusters = getLeafObjsIDs(node, Clusters);
k = 0;
for c = 1: length(Clusters)
    k = k + 1;
    pi(Clusters{c},1) = k;
end
end

%%
function Clusters = getLeafObjsIDs(node, Clusters)
% Base case: if node is empty, return
if isempty(node)
    return;
end
% If it is a leaf node, add its objs_ID to the cell array
if ~isstruct(node.left) && ~isstruct(node.right)
    Clusters{end + 1, 1} = node.objs;
    return;
end
% Otherwise, continue the recursion
if isstruct(node.left)
    Clusters = getLeafObjsIDs(node.left, Clusters);
end
if isstruct(node.right)
    Clusters = getLeafObjsIDs(node.right, Clusters);
end
end