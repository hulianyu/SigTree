% Draw the tree
function drawTree(node, rowNames, I)
% Calculate the tree's depth and the number of leaf nodes
maxDepth = treeDepth(node);
numLeaves = countLeaves(node);
% Set initial parameters for drawing the tree
initialX = 0.5;
initialY = 1;
deltaX = 1 / (2 ^ (maxDepth - 1));
deltaY = 0.8 / maxDepth;
% Draw the tree
figure;
hold on;
axis off;
title(sprintf('%s - Max Depth: %d, Leaves: %d', rowNames{I}, maxDepth, numLeaves)); % Set the title
drawNode(node, initialX, initialY, deltaX, deltaY);
hold off;
end

% Count the number of leaf nodes
function leafCount = countLeaves(node)
if isempty(node)
    leafCount = 0;
elseif ~isstruct(node.left) && ~isstruct(node.right)
    leafCount = 1;
else
    leafCountLeft = countLeaves(node.left);
    leafCountRight = countLeaves(node.right);
    leafCount = leafCountLeft + leafCountRight;
end
end

% Draw each node of the tree
function drawNode(node, x, y, deltaX, deltaY)
if isempty(node)
    return;
end
% Check if node.category has sufficient elements
if length(node.category) >= 2
    nodeText = sprintf('p:%.0e\nm:%d\nq:%d\nN:%d', ...
        node.pval, node.category(1), node.category(2), length(node.objs));
else
    % Use default values or alternative logic if node.category does not have two elements
    nodeText = sprintf('N:%d', ...
        length(node.objs));
end
text(x, y, nodeText, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
if ~isempty(node.left) && isstruct(node.left)
    newX = x - deltaX;
    newY = y - deltaY;
    line([x, newX], [y, newY], 'Color', 'black');
    drawNode(node.left, newX, newY, deltaX / 2, deltaY);
end
if ~isempty(node.right) && isstruct(node.right)
    newX = x + deltaX;
    newY = y - deltaY;
    line([x, newX], [y, newY], 'Color', 'black');
    drawNode(node.right, newX, newY, deltaX / 2, deltaY);
end
end

% % Calculate the depth of the tree
% function maxDepth = treeDepth(node)
% if isempty(node) || (~isstruct(node.left) && ~isstruct(node.right))
%     maxDepth = 0;
% else
%     maxDepthLeft = treeDepth(node.left);
%     maxDepthRight = treeDepth(node.right);
%     maxDepth = 1 + max(maxDepthLeft, maxDepthRight);
% end
% end

% function node = createNode(objs_ID, min_pval, best_mcat, left, right)
% % Create a new tree node containing min_pval, left child, and right child
% node.objs = objs_ID;
% node.pval = min_pval;
% node.category = best_mcat;
% node.left = left;
% node.right = right;
% end