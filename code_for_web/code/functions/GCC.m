function gcc = GCC(A)
%ONLYGCC Summary of this function goes here
%   Detailed explanation goes here
[cb, sizes] = conncomp(graph(A));
[~,indGcc] = max(sizes);
gcc = find(cb==indGcc);
end

