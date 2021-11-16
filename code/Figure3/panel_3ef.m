%% Maps of approximations for the shells mean field

addpath('..\functions');

%% map - 1st 2nd and 3rd order lower bound
% global a b

NameOfModel = 'MM'; % 'Neural'; % 'Glauber'; %  'Eco'; % 'Simple'; %  'Voter'; % 'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );

% a=1; b=2; % exponents of MM

w = 0.7;
k = 10;

colors = [ 255 0 0; 0 210 200]/255;

approx_ord = 1; % 
A = -0.01; B = 10;
x = setdiff(linspace(A,B,201),0);

% find the states of a free system using mean field
f_MF = @(x) M{1}(x)+w*(k+1)*M{2}(x).*M{3}(x); 
fpMF = arrayfun(@(x0) fzero(f_MF,x0),x);
fpMF = uniquetol(fpMF,1e-2,'DataScale',1);
x_low = fpMF(1); x_high = fpMF(3);

A=x_low-1/2; B=x_high+1; % limits of the interval
x=linspace(A,B,201);

R = @(x) M{1}(x)./M{2}(x); 
M2 = M{3};
% Rinv = find_inv_func_for_nnmono(@(x) -F(x),setdiff( linspace(A,B,1001),0));
Rinv = M{5};
switch approx_ord
    case 0        
        F = @(x) -R(x)/w-(k)*M2(x_low);
    case 1
        F = @(x) -R(x)/w-(k)*M2(Rinv(w*M2(x)+w*(k)*M2(x_low)));
end
[F_mon, xmax , ymax] = find_upbnd_func_for_nnmono2(F,setdiff( linspace(A,B,1001),0));

figure; hold on;

plot(x,arrayfun( F_mon,x),'color',[153, 0, 204]/256,'LineWidth',3)
plot(x,F(x),'--','color',[153, 0, 204]/256,'LineWidth',1.5)
plot(x,M2(x),'color',[255, 153, 0]/256,'LineWidth',3)

ax = gca;
ax.XTick = [];
xlabel('$x$','Interpreter','latex')
set(gca,'YTick',[],'FontSize',30,'Box','on','LineWidth',1.5)
set(gcf,'Color','w')
axis([x(1),2.5,-0.7,2])
xlim([-0.3,2.8])
ylim([-0.1,1.1])

f_s = @(x) F_mon(x)-M2(x);
% fp1 = arrayfun(@(x0) fsolve(f_s,x0),x);
fp = arrayfun(@(x0) fzero(f_s,x0),x);
% fp = union(fp1,fp2);
fp(isnan(fp))=[];
fp = uniquetol(fp,1e-2,'DataScale',1);

if  length(fp)>2
    plot([fp(2),fp(2)],ylim,'--k','LineWidth',1.5)
    ax.XTick = fp(2);
    ax.XTickLabel = {'$\Delta_c$'}; 
    scatter(fp([1,3]),M2(fp([1 3])),90,colors,'filled')
    scatter(fp(2),M2(fp(2)),90,'k','MarkerFaceColor','w')
elseif length(fp)==1
    scatter(fp,M2(fp),90,'r','filled')
end
ax.XTick = [];


%% propagation
steps = 5;

j=0;
for x0 = fp(2)+0.15*[-1,1]
    j=j+1;
    xl = zeros(1,1+steps);
    xl(1) = x0; %+bias	
    for i=1:steps
        xl(i+1)= fzero(@(x) F_mon(x)-M2(xl(i)),[x_low,x_high]);
    end
    xlp = repelem(xl,2);
    yl = M2(xlp);
    yl(1:2:end) = F(xlp(1:2:end));
    yl(1) = ax.YLim(1);

    for i=1:2:length(xlp)
        hh = plot(xlp(1:i),yl(1:i),'LineWidth',2,'Color',colors(j,:));       
    end
    
end


%% to be on the top

if  length(fp)>2
     scatter(fp([1,3]),M2(fp([1 3])),90,colors,'filled')
    scatter(fp(2),M2(fp(2)),90,'k','MarkerFaceColor','w')
elseif length(fp)==1
    scatter(fp,M2(fp),90,'r','filled')
end

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'ef';
save_pdf_min_size([folder,filename])

