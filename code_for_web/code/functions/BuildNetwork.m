function A = BuildNetwork( n, NetStruct , parameters ,only_gcc)
%BUILDNETWORK builds a network according to input parameters
%   options - 'ER','SF','Lattice', 'RR', 'RT', 'complete'.
%   parameters = k (in 'ER','RR') / [lambda ,k0] (at 'SF')

switch NetStruct
    case 'ER'
        for ER=0 % Erdos-Renyi
            k = parameters;
            p = k/(n-1);
            A = logical(sprandsym(n,p));
            A(1:n+1:end) = 0;
        end
    case 'SF'
        gamma = parameters(1);
        k0 = parameters(2);
        A = BuildSF(n,gamma,k0);
    case 'Lattice'
        for Lattice=0 % Lattice
            L = floor(sqrt(n));
            I = speye(L);
            E = sparse(2:L,1:L-1,1,L,L);
            if strcmp(parameters,'periodic'); E(L)=1; end
            D = E + E';
            A = kron(I,D) + kron(D,I);
        end
    case 'RR' % Random-Regular
        k = parameters;
        A = BuildRR( n,k );
    case 'RRpure' % pure Random-Regular
        k = parameters;
        A = BuildRRpure( n,k );
    case 'RRMultDeg' % Random-Regular with 2 degrees
        ks = parameters;        
        A = BuildRRMultDeg( n,ks );        
    case 'RT' % Regular-Tree
        k = parameters;
        A = BuildRegTree( n,k );
    case 'complete'
        A = ones(n,n);
        A(1:n+1:end)=0;
    case 'CM'
        P = parameters;
        A = BuildCM(n,P);
    case {'PPI_Human' , 'PPI_Yeast' , 'Brain', 'Brain_perc' , 'ECO1', 'ECO1_perc', 'ECO2', 'ECO2_perc'}
        folder = '..\..\data\';
        st = load([folder,NetStruct,'.mat']);
        A = st.A;
end

% remove the isolated nodes outside the gcc
if nargin == 4
    if strcmp( only_gcc , 'gcc')
        A = onlyGCC(A);
    else
        disp('warning: unknown input, should be ''gcc'' in order to remove isolated nodes')
    end
end

