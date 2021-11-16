function set_weights_to_Aij(DistName, parameters)
global A Anw

A = double(Anw);
m = nnz(Anw);

switch DistName
    case 'normal'
        mu = parameters(1);
        sig = parameters(2);
        w = mu+sig*randn(m,1);
        w(w<0) = 0;
    otherwise
        disp(['what is ',DistName,'??'])
end

A(Anw)=w;
