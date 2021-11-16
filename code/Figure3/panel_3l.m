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

N=1e4;
wVec = logspace(-1,log10(20),50);
k0 = 2;
DeltaVec = logspace(-2.5,0.5,50);

all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'

kappaEdge = 15+1/2*[-1,1];

reals = 20;

phaseMat = zeros(length(DeltaVec),length(wVec));

Xss = zeros(1,2);

for id = 1:length(DeltaVec)
    Delta = DeltaVec(id)
    
    for iw=1:length(wVec)
        w = wVec(iw);
        
        ik = 0;
        while ik<reals
            
            NetStruct = 'SF';
            c = rand<0.3;
            gamma = c*(rand+3) + (1-c)*(0.5*rand+2.5);
            parameters = [gamma k0];
            
            Anw = BuildNetwork(N, NetStruct, parameters,'gcc'); % Adjacency matrix, not weighted
            n = size(Anw,1);
            
            k = mean(sum(Anw));
            k2 = mean(sum(Anw).^2);
            kappa = k2/k-1;
            
            isin = discretize(kappa,kappaEdge);
            if isnan(isin); continue; end
            ik = ik+1;
            
            A = w*Anw;
            Factor_Xeff = sum(A,1)/sum(sum(A));
            
            % initial/fixed condition
            x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
            toplot = 0;
            
            for idx_init=1:length(free_value_vec)
                conditions.free_value = free_value_vec(idx_init);
                if idx_init==2 && w>0.4
                    Xss(idx_init) = 10;
                else
                    Xss(idx_init) = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release);
                end
            end
            phase = sum(Xss>x_th);
            phaseMat(id,iw) = phaseMat(id,iw) + phase;
        end
    end
end
phaseMat = phaseMat./reals;

toc

%% Theory
disp('theory')
options = optimset('TolX',1e-3);
dc_theory = arrayfun(@(w) find_dc_by_wk_theory_intersection(w,mean(kappaEdge),range) ,  wVec);

toc

%% figure

figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'XScale','log','YScale','log','layer','top')
xlabel('\boldmath$\omega$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\Delta$','Interpreter','latex','FontSize',30)
axis tight square
xticks(10.^(-2:2))
yticks(10.^(-2:2))

imagesc(wVec,DeltaVec,phaseMat)

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
colorbar('LineWidth',2,'Ticks',[0,1,2])

plot(wVec,dc_theory,'w','LineWidth',2)
idxc_unrec_theo = find(isfinite(dc_theory),1);
ylims = ylim; ymax = ylims(2);
plot(wVec([idxc_unrec_theo, idxc_unrec_theo]),[dc_theory(idxc_unrec_theo),ymax ],'w','LineWidth',2);

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'l';
save_pdf_min_size([folder,filename])


