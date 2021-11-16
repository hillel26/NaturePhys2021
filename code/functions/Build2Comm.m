function [A] = Build2Comm(A1,A2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
n1 = size(A1,1);
n2 = size(A2,1);
A = blkdiag(A1,A2);
s = randi(n1);
t = n1 + randi(n2);
A(s,t) = 1;
A(t,s) = 1;
end

