function Xss = find_x_via_omega(omega,M,all_or_eff,toplot,conditions,release)
global A Anw x0 Factor_Xeff
A = omega*Anw;
Factor_Xeff = sum(A,1)/sum(sum(A));
Xss = SolveOdes(x0,M,all_or_eff, toplot,conditions,release);