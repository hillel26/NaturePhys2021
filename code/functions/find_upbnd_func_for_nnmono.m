function [fmon,x,y] = find_upbnd_func_for_nnmono(f,x)
% consider a nonmonotonic function,
% find upper monotonic increasing bound
% inputs 
% f - input function
% x - vector of x values in the region where the fmon should be calculated

y = 0*x;
y(1) = f(x(1));
for i=2:length(x)
    y(i) = max(y(i-1),f(x(i)));    
end
fmon = @(xq) interp1(x,y,xq);