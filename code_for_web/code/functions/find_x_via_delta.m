function Xss = find_x_via_delta(delta,M,all_or_eff,toplot,conditions,release)
global x0
x0 = delta*logical(x0);
Xss = SolveOdes(x0,M,all_or_eff, toplot,conditions,release);