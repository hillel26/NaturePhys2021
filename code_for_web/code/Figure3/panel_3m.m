clear; format compact; clc;

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
free_value_vec = [1e-3,10]; % the initial value of the free nodes
NumExcited = 1;
release = 1; % to release or not after holding

N = 1e4;
w = 0.7;
DeltaVec = logspace(-1,1,50);

all_or_eff = 'eff'; % which xss we need: 'all'/ 'eff' /'eff_F'

kappaEdges = linspace(1,50,50);
kappaCounter = zeros(1,length(kappaEdges)-1);
kappaVec = ( kappaEdges(1:end-1)+kappaEdges(2:end))/2;

many_reals = 20; % number of realizations
few_reals = 20;
reals = ((DeltaVec>0.3)&(DeltaVec<2))' | ((kappaVec>0)&(kappaVec<15));
reals = reals*(many_reals-few_reals)+few_reals;

phaseMat = zeros(length(DeltaVec),length(kappaVec));

Xss = zeros(1,2);

for id = 1:length(DeltaVec)
    Delta = DeltaVec(id)
        
    kappaCounter = 0*kappaCounter;
    while any(kappaCounter<reals(id,:)) 
                
        if rand<0.95
            NetStruct = 'SF';
            randNum = rand<0.3;
            gamma = randNum*(rand+3) + (1-randNum)*(0.5*rand+2.5);
            k0 = rand*1.5+1.5;
            parameters = [gamma k0];
        else
            NetStruct = 'ER';
            k0 = rand*2+1;
            parameters = k0;
        end
        
        Anw = BuildNetwork(N, NetStruct, parameters,'gcc'); % Adjacency matrix, not weighted
        n = size(Anw,1);
        
        k = mean(sum(Anw));
        k2 = mean(sum(Anw).^2);
        kappa = k2/k-1;
        
        ik = discretize(kappa,kappaEdges);
        if isnan(ik) || kappaCounter(ik) == reals(id,ik); continue; end
        kappaCounter(ik) = kappaCounter(ik)+1;
        
        A = w*Anw;
        Factor_Xeff = sum(A,1)/sum(sum(A));
        
        % initial/fixed condition
        x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
        toplot = 0;
        
        for idx_init=1:length(free_value_vec)
            conditions.free_value = free_value_vec(idx_init);
            if idx_init==2 && kappa>5
                Xss(idx_init) = 10;
            else
                Xss(idx_init) = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
            end
        end
        phase = sum(Xss>x_th);
        phaseMat(id,ik) = phaseMat(id,ik) + phase;
    end
end
phaseMat = phaseMat./reals;

toc

%% Theory
disp('theory')
options = optimset('TolX',1e-3);
dc_theory = arrayfun(@(kappa) find_dc_by_wk_theory_intersection(w,kappa,range) , kappaEdges);

toc

%% figure

figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'XScale','lin','YScale','log','layer','top')
xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\Delta$','Interpreter','latex','FontSize',30)
axis tight square
xticks([10 40])
yticks([1e-1 1e0 1e1])

imagesc(kappaVec,DeltaVec,phaseMat)

clrs = [0 51 102; 255 192 0; 155 0 0]/255;
ia_clr =  [155 0 0]/255;
urec_clr = [255 192 0]/255;
rec_clr = [0 51 102]/255;

alpha = linspace(0,1,10)';
map1 = (1-alpha)*ia_clr + (alpha)*urec_clr;
map2 = (1-alpha)*urec_clr + (alpha)*rec_clr;
map = [map1;map2];
colormap( map)
caxis([0,2])
% colorbar('LineWidth',2,'Ticks',[0,1,2])

plot(kappaEdges,dc_theory,'w','LineWidth',2)
idxc_unrec_theo = find(isfinite(dc_theory),1);
ylims = ylim; ymax = ylims(2);
plot(kappaEdges([idxc_unrec_theo, idxc_unrec_theo]),[dc_theory(idxc_unrec_theo),ymax ],'w','LineWidth',2);

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'm';
save_pdf_min_size([folder,filename])


