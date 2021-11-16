function [A] = BuildCM(N,P)
% generate general network by the configuration model given degree distribution
% P is the function of the degree distribution p(k)
%%

nk = P(1:N)/sum(P(1:N))*N;% # of nodes with degree k 

for i=N:-1:2    
    nk(i-1) = nk(i-1) + mod(nk(i),1);
    nk(i) = floor(nk(i));
end
nk(1) = round(nk(1));

ki = repelem(1:N,nk);

stubs = repelem(1:N,ki);
stubs = stubs(randperm(length(stubs)));

if mod(length(stubs),2) % if it's odd
    stubs(end)=[];
end

s = stubs(1:2:length(stubs));
t = stubs(2:2:length(stubs));
A = sparse([s t],[t s],1,N,N); % Symetric adjacency matrix
A = logical(A); % Because multiplicity of links. yeah..
A(1:size(A,1)+1:end) = 0;


