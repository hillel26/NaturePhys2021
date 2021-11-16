clear; format compact;% close all; clc;

addpath('..\functions');

global A Anw Factor_Xeff NameOfModel

tic

% model of dynamics
NameOfModel = 'MM'; % 'Neural'; %  'Eco'; % 'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );
range = [0,100]; % range of values of x for theory
x_th = 1e-2;

conditions.type = 'BC';
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
Delta = 10;
free_value_vec = [1e-3,10]; % the initial value of the free nodes
NumExcited = 1;
release = 1; % to release or not after holding

N = 1e4;
% k0 = 2;
wVec = logspace(-1.5,0.5,50); % linspace(0.1,1,50); % 

all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'

many_reals = 20; % number of realizations
few_reals = 1;

kappaEdges = linspace(1,50,50); % linspace(5,10,50); % 
kappaCounter = zeros(1,length(kappaEdges)-1);
kappaVec = ( kappaEdges(1:end-1)+kappaEdges(2:end))/2;

phaseMat = zeros(length(wVec),length(kappaVec));

Xss = zeros(1,2);

for iw = 1:length(wVec)
    w = wVec(iw)
    
    if w<1.5 && w>0.03; reals = many_reals; else; reals = few_reals; end % in order to not waste time on trivial areas
    
    kappaCounter = 0*kappaCounter;
    while any(kappaCounter<reals)
        
        NetStruct = 'SF';
        a = rand<0.3;
        gamma = a*(rand+3) + (1-a)*(0.5*rand+2.5);
        k0 = rand*2+1;
        parameters = [gamma k0];
        
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
        
        % initial/fixed condition
        x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
        toplot = 0;
        
        for idx_init=1:length(free_value_vec)
            conditions.free_value = free_value_vec(idx_init);
            
            Xss(idx_init) = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
        end
        phase = sum(Xss>x_th);
        phaseMat(iw,ik) = phaseMat(iw,ik) + phase;
    end
    phaseMat(iw,:) = phaseMat(iw,:)/reals;
end

toc

%% Theory
disp('theory')
options = optimset('TolX',1e-3);
wc_theory = arrayfun(@(kappa) fzero(@(w) is_rcoverable_theory(w,kappa,Delta,range) - 1/2 , [1e-1,1e1],options) , kappaEdges);

toc

%% figure

figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'XScale','lin','YScale','lin','layer','top')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\omega$','Interpreter','latex','FontSize',30)
axis tight square
xticks([5, 10])
yticks([0.5,1])

imagesc(kappaVec,wVec,phaseMat)

clrs = [0 51 102; 255 192 0; 155 0 0]/255;
ia_clr =  [155 0 0]/255;
urec_clr = [255 192 0]/255;
rec_clr = [0 51 102]/255;

alpha = linspace(0,1,30)';
map1 = (1-alpha)*ia_clr + (alpha)*urec_clr;
map2 = (1-alpha)*urec_clr + (alpha)*rec_clr;
map = [map1;map2];
colormap( map)
caxis([0,2])
% colorbar('LineWidth',2,'Ticks',[0,1,2])

plot(kappaEdges,wc_theory,'w','LineWidth',2)


%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'k';
save_pdf_min_size([folder,filename])
