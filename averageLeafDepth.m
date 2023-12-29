% Calculate the average depth of leaf nodes
function avgLeafDepth = AverageLeafDepth(node)
    % Call the helper function to get the sum of depths and the number of leaves
    [sumOfDepths, numLeaves] = leafDepthsHelper(node, 0);

    % Calculate the average leaf depth
    avgLeafDepth = sumOfDepths / numLeaves;
end

% Helper function to calculate the sum of depths of leaf nodes
function [sumOfDepths, numLeaves] = leafDepthsHelper(node, currentDepth)
    % Initialize sumOfDepths and numLeaves
    sumOfDepths = 0;
    numLeaves = 0;

    % Check if the current node is empty
    if isempty(node)
        return;
    end

    % If the current node is a leaf node
    if ~isstruct(node.left) && ~isstruct(node.right)
        sumOfDepths = currentDepth;
        numLeaves = 1;
    else
        % Traverse left and right subtrees
        [sumLeft, numLeft] = leafDepthsHelper(node.left, currentDepth + 1);
        [sumRight, numRight] = leafDepthsHelper(node.right, currentDepth + 1);

        % Accumulate the sums and counts
        sumOfDepths = sumLeft + sumRight;
        numLeaves = numLeft + numRight;
    end
end