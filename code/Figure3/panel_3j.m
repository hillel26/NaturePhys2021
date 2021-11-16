clc; clear; format compact;
% close all; 

addpath('..\functions');

global NameOfModel

tic

% model of dynamics
NameOfModel = 'MM'; % 'Eco'; % 'Simple'; % 'Glauber'; % 'Voter'; % 'SIS';% 'MAK';%  'PD';%
range = [0,100];

wVec = linspace(0.1,2,50);
kappaVec = linspace(2,50,100);

%% theory

dc_theory = zeros(length(wVec),length(kappaVec));

for iw = 1:length(wVec)
    w = wVec(iw)
    
    for ik = 1:length(kappaVec)
        if ~mod(ik,10); fprintf('%d ',ik); end
        kappa = kappaVec(ik);
        
        dc_theory(iw,ik) = find_dc_by_wk_theory_intersection(w,kappa,range);
    end
    fprintf('\n')
end


%% figure

figure; hold on

set(gca,'FontSize',20,'box','on','LineWidth',2,'layer','top','XScale','lin','YScale','lin')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\omega$','Interpreter','latex','FontSize',30)
zlabel('\boldmath$\Delta$','Interpreter','latex','FontSize',30)
axis square
axis([0 50 0 2 0 2.5])
xticks(0:20:50); yticks(0:2); zticks(0:3)

K = repmat(kappaVec,length(wVec),1);
W = repmat(wVec',1,length(kappaVec));

% dc_theory(dc_theory==0) = inf;
dc_theory0 = dc_theory;
dc_theory0(isinf(dc_theory0)) = 0;

surf(kappaVec,wVec,dc_theory0,'EdgeColor','none','EdgeAlpha',0.1,'FaceAlpha',0.5) % ,'EdgeAlpha',0.3
% patch()
map = [0.6 0 0; parula];
colormap(map)

% equal-curves
ws = 0.7;
[~, iws] = min(abs(wVec-ws));
ks = 15;
[~, iks] = min(abs(kappaVec-ks));
% ds = 1;

plot3(K(iws,:),W(iws,:),dc_theory(iws,:),'color',[150 0 230]/255,'LineWidth',3)
plot3(K(:,iks),W(:,iks),dc_theory(:,iks),'color',[0 220 50]/255,'LineWidth',3) 
% contour3(K,W,dc_theory,[ds ds],'LineColor','k','LineWidth',3)

wc_inds = arrayfun ( @(i) find(isinf(dc_theory(:,i)),1,'last'), 1:length(kappaVec) );
wc = wVec(wc_inds);
plot(kappaVec,wc,'color',[0 0 250]/255,'LineWidth',3)

view(230,25)

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'j';
save_pdf_min_size([folder,filename])

