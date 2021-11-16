function dc = find_dc_by_wk_theory_directed_intersection(w,kappa,range)
% 2nd order lower bound
global  NameOfModel

% inf_val = 1e3; % for helping search with fzero
% w=3; k=5;
bnd_ord = 0;
A = range(1); B = range(2);
x = A+logspace(-3,log10(B),1e3);

M = KindOfDynamics( NameOfModel );

% find the states of free system
f_MF = @(x) M{1}(x)+w*(kappa+1).*M{2}(x).*M{3}(x);
if ~f_MF(A); x = [A,x]; end % to include 0 for MM
options = optimset('display','off','TolX',1e-10);
fpMF = arrayfun(@(x0) fzero(f_MF,x0,options),x);
fpMF = unique(round(fpMF,5));

[~, xmax] = findpeaks(f_MF(x),x);

if length(fpMF)<2
    if strcmp(NameOfModel,'Eco')
         dc = 0;
        return
    elseif ( fpMF(1) > xmax )
        dc = 0;
        return
    else
        dc = inf;
        return
    end
end

x_low = fpMF(1); % x_high = fpMF(3);
% A = x_low-1/2; B = x_high; % limits of the interval

R = @(x) -M{1}(x)./M{2}(x);
M2 = M{3};
rhs = @(x) M2(x);
if length(M) > 4
    Rinv = M{5};
else
    Rinv = find_inv_func_for_nnmono(@(x) R(x),x);
end

switch bnd_ord
    case 0
        F = @(x) R(x)/w-(kappa)*M2(x_low);
    case 1
        F = @(x) R(x)/w-(kappa)*M2(Rinv(w*M2(x)+w*(kappa)*M2(x_low)));
    case 2
        xl2 = @(xl,xl2) w*M2(Rinv(w*M2(xl)+w*(kappa)*M2(xl2)));
        F = @(x) R(x)/w-1/w*(kappa)*xl2(x,xl2(x,xl2(x,xl2(x,xl2(x,xl2(x,x_low))))));
    case 3
        g1 = @(x) Rinv(w*M2(x)); g2 = @(x) Rinv(w*(kappa)*M2(x));
        rhs = @(x) g2(g1(g1(x)));
        F = @(x) x;
end
F_mon = find_upbnd_func_for_nnmono(F,x);

f_s = @(x) F_mon(x)-rhs(x);
% f_s = @(x) F(x)-rhs(x); %  caution !!! test. the line above should be true  
fp = arrayfun(@(x0) fzero(f_s,x0,options),x);
fp(isnan(fp))=[];
fp = unique(round(fp,5));

if length(fp)>1
    dc = fp(2);
else
    dc = inf;
end

% figures of functions
% plot(x,arrayfun(F_mon,x),'--k','LineWidth',1.5)
% plot(x,F(x),'color',[153, 0, 204]/256,'LineWidth',2)
% plot(x,rhs(x),'color',[255, 153, 0]/256,'LineWidth',2)
%
% ax = gca;
% ax.XTick = [];
% xlabel('$x$','Interpreter','latex')
% set(gca,'YTick',[],'FontSize',20,'Box','on','TickLabelInterpreter','latex')
% set(gcf,'Color','w')
% axis tight


% if  length(fp)>2
%     plot([fp(2),fp(2)],ylim,'--k','LineWidth',1)
%     ax.XTick = fp(2);
%     ax.XTickLabel = {'$\Delta_c$'};
%     scatter(fp([1,3]),rhs(fp([1 3])),50,'k','filled')
%     scatter(fp(2),rhs(fp(2)),50,'k','MarkerFaceColor','w')
% elseif length(fp)==1
%     scatter(fp,rhs(fp),50,'k','filled')
% end
% ax.XTick = [];
% text(fp(1),lhs_mon(fp(2)),'$F(x)$','Interpreter','latex','FontSize',15,'VerticalAlignment','bottom')
% text(1.3,rhs(1.3),'$M_2(x)$','Interpreter','latex','FontSize',15,'VerticalAlignment','top')
