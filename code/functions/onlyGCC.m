function [A, gcc_inds] = onlyGCC(A)
%ONLYGCC Summary of this function goes here
%   Detailed explanation goes here
[cb, sizes] = conncomp(graph(A));
[~,indGcc] = max(sizes);
gcc = find(cb==indGcc);
A = A(gcc,gcc);

if nargout>1
   gcc_inds =  gcc;
end
end

