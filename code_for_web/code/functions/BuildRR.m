function [ A ] = BuildRR( N,k )

%BUILDRR generate random regular graph
stubs = repelem(1:N,k);

stubs = stubs(randperm(length(stubs)));

s = stubs(1:2:end);
t = stubs(2:2:end);

Ain = sparse(s,t,1,N,N); % Adjacency matrix
A = logical(Ain +Ain');

A(1:size(A,1)+1:end) = 0;
A = A(sum(A,2)>0,:);
A = A(:,sum(A,1)>0);


