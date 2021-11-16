%% generate SF (or other degree distribution) network by the configuration model
function [A] = BuildSF(N,lambda,k0,kf)


%%
% kf = k0*N^(1/(lambda-1));

if nargin>3
    ki = round( k0 * (1 + ((kf/k0).^(1-lambda)-1)*rand(1,N)).^(1/(1-lambda)) );
else
    ki = round( k0 * rand(1,N).^(1/(1-lambda)) );
end

stubs = repelem(1:N,ki);
stubs = stubs(randperm(length(stubs)));
if mod(length(stubs),2) % if it's odd
    stubs(end)=[];
end

I = stubs(1:2:length(stubs));
J = stubs(2:2:length(stubs));
A = sparse([I J],[J I],1,N,N); % Symetric adjacency matrix
A = logical(A); % Because multiplicity of links. yeah..
A(1:size(A,1)+1:end) = 0;

%% wrong way!!!!
% k0 = round(k0);
% P = @(k) k.^(-lambda); 
% kf = N;
% nk = zeros(1,N);
% nk(k0:kf) = P(k0:kf)/sum(P(k0:kf))*N;% # of nodes with degree k 
% for i=kf:-1:k0+1    
%     nk(i-1) = nk(i-1)+mod(nk(i),1);
%     nk(i) = floor(nk(i));
% end
% kf = find(nk, 1, 'last' ); % because of round kf changes
% ki = repelem(1:n,nk);

