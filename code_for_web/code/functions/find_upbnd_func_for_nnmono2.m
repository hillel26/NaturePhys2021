function [fmon,xmax,ymax] = find_upbnd_func_for_nnmono2(f,x)
% consider a nonmonotonic function,
% find upper monotonic increasing bound
% inputs
% f - input function
% x - vector of x values in the region where the fmon should be calculated

[ymax, xmax] = findpeaks(f(x),x);

% [xmax,ymax] = fminsearch(@(x) -f(x),0.2);

if isempty(xmax)
    fmon = f;
else
    fmon = @fmono;
end


    function out = fmono(x)
        out = f(x);
        if x>xmax
            out = max(ymax, f(x));
        end
    end


end