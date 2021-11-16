function [xss,tau1,tau2] = find_x_and_tau_via_omega(omega,x0,M,toplot,conditions_type,times)
global A Anw Factor_Xeff
A = omega*Anw;
Factor_Xeff = sum(A,1)/sum(sum(A));
[xss,tau1,tau2] = SolveOdes(x0,M,'eff', toplot,conditions_type,times);