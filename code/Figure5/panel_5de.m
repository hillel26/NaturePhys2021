%% introductions

clear; clc; format compact;% close all; 

addpath('..\functions');

global A Anw Factor_Xeff NameOfModel mu

tic

%% model of dynamics
NameOfModel = 'Neural'; %  'MM'; % 'Eco'; % 'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
mu = 10; % parameter of brain dynamics
M = KindOfDynamics( NameOfModel );
range_x = [0,300]; % range of values of x for theory
x_th = 1; % threshold

%% reviving details
conditions.type = 'BC';
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
Delta = 8;
free_value_high = 10;
free_value_low = 0;

release = 1; % to release or not after holding

%% simulations

N=1e4;
wVec = logspace(log10(0.4),3,50);

all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'

reals = 10; % number of realizations

kappaEdges = linspace(2,30,51);
kappaCounter = zeros(1,length(kappaEdges)-1);
kappaVec = ( kappaEdges(1:end-1)+kappaEdges(2:end))/2;

phaseMat = zeros(length(wVec),length(kappaVec));

for iw = 1:length(wVec)
    w = wVec(iw)
        
    kappaCounter = 0*kappaCounter;
    while any(kappaCounter<reals)
        
        if rand<0.9
            NetStruct = 'SF';
            randNum = rand<0.3;
            gamma = randNum*(rand+3) + (1-randNum)*(0.5*rand+2.5);
            k0 = rand*1.5+1.5;
            parameters = [gamma k0];
        else
            NetStruct = 'ER';
            k0 = rand*3+1;
            parameters = k0;
        end        
        
        Anw = BuildNetwork(N, NetStruct, parameters,'gcc'); % Adjacency matrix, not weighted
        n = size(Anw,1);
        
        k = mean(sum(Anw));
        k2 = mean(sum(Anw).^2);
        kappa = k2/k-1;
        
        ik = discretize(kappa,kappaEdges);
        if isnan(ik) || kappaCounter(ik) == reals; continue; end
        kappaCounter(ik) = kappaCounter(ik)+1;
        
        A = w*Anw;
        Factor_Xeff = sum(A,1)/sum(sum(A));
        
        toplot = 0;
        
        % high free
        NumExcited = 0;
        conditions.free_value = free_value_high;
        x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
        Xss = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
        if Xss<x_th
            phase = 0; % inactive
        else
            % low free
            NumExcited = 0;
            conditions.free_value = free_value_low;
            x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
            Xss = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
            if Xss>x_th
                phase = 3; % active
            else
                % ignition
                NumExcited = 1;
                conditions.free_value = free_value_low;
                x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
                Xss = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
                if Xss>x_th
                    phase = 2; % recoverable
                else
                    phase = 1; % unrecoverable 
                end
            end
        end
        phaseMat(iw,ik) = phaseMat(iw,ik) + phase;
    end
    phaseMat(iw,:) = phaseMat(iw,:)/reals;
end

toc

%% Theory
tic
disp('theory')
kappa_theory = kappaEdges; % linspace(1,30,40);
wc_theory = 0*kappa_theory;

options = optimset('TolX',1e-4,'Display','off');
for ik=1:length(kappa_theory)
    kappa = kappa_theory(ik)
    wc_theory(ik) = fzero(@(w) is_rcoverable_theory(w,kappa,Delta,range_x) - 1/2 , ...
        [1e-1,1e2],options) ;
end

toc

%% figure

figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'XScale','lin','YScale','log','layer','top')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\omega$','Interpreter','latex','FontSize',30)
axis tight square
xticks(0:10:40); yticks(10.^(-2:2:3))

imagesc(kappaVec,wVec,phaseMat)

ia_clr =  [155 0 0]/255;
urec_clr = [255 192 0]/255;
rec_clr = [0 51 102]/255;
ac_clr = [0 150 150]/255;

% urec_clr = [200 200 200]/255;
% rec_clr = [200 200 200]/255;

alpha = linspace(0,1,10)';
map1 = (1-alpha)*ia_clr + (alpha)*urec_clr;
map2 = (1-alpha)*urec_clr + (alpha)*rec_clr;
map3 = (1-alpha)*rec_clr + (alpha)*ac_clr;

map = [map1;map2;map3];
colormap( map)
caxis([0,3])
% colorbar('LineWidth',2,'Ticks',[0,1,2])

axiss = [xlim,ylim];
plot(kappa_theory,wc_theory,'w','LineWidth',2)
axis(axiss)


%% save the figure in the folder 'output'
folder = '..\..\output\Figure5\';
filename = 'de';
save_pdf_min_size([folder,filename])

