function [finv,x_unq,y_unq] = find_inv_func_for_nnmono(f,x)
% consider a nonmonotonic function,
% find upper monotonic increasing bound
% find its inverse function
% inputs 
% f - input function
% x - vector of x values in the region where the finv should be calculated

y = 0*x;
y(1) = f(x(1));
for i=2:length(x)
    y(i) = max(y(i-1),f(x(i)));    
end
[y_unq, ind] = unique(y);
x_unq = x(ind);
finv = @(x) interp1(y_unq,x_unq,x);