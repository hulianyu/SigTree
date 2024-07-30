function n = simplifiedHistcounts(x, edges)
% Simplified histogram bin counts for pre-defined edges.
% This version assumes the input x is valid, non-NaN, non-Inf double precision numbers,
% and the edges are pre-defined and valid.
% Calculate histogram counts directly with pre-defined edges
n = matlab.internal.math.histcounts(x, edges);
end
