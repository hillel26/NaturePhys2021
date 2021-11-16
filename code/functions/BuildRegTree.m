function [ A ] = BuildRegTree( n,k )

% generate k-regular tree with n nodes around the source 1
kiout = repelem([k,k-1],[1,n-1]);
s = repelem(1:n, kiout);
t = 2:n;

Aout = sparse(s(1:n-1),t,1,n,n); % Adjacency matrix
A = logical(Aout +Aout');
