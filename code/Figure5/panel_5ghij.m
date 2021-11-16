%% w_inter vs w_intra diagram
%% introductions

clear; format compact;

global A Anw x0 A1 A2 A12 A21 conditions M comm1 mu

addpath('..\functions');

load('..\..\data\Brain')

tic

%% dynamics
NameOfModel = 'Neural'; % 'Eco'; % 'MM'; %  'Simple'; % 'Glauber'; % 'Voter'; % 'SIS';% 'MAK';%  'PD';%
mu = 10;

h = gobjects(0);
M = KindOfDynamics( NameOfModel );
range_sols = [0,50]; % range of high and low state for the search of the theory solution

ChooseExcNodes = 'specific'; % options - {'low degree', 'high degree', 'local', 'rand','1'}
conditions.type = 'BC';% 'IC'; %
release = 0;
free_value_vec = [0, 10];
all_or_eff = 'average'; % which xss we need: 'all'/ 'eff' /'eff_F'/ 'average'

diff_conds = {'ignited'}; %{'free low','free high','ignited'};
th1 = 10;
th2 = 20;

Delta = 20; % logspace(-1,2,4);
NExcVec = [0 1]; %unique( round(logspace(0,log10(NVec/25),21)) );

reals = 20;
win_c_reals1 = zeros(1,reals);
win_c_reals2 = win_c_reals1;
xReals = zeros(1,reals);
stateReals = zeros(1,reals);

%% build the network
Anw = (A>0.03);
Anw = sparse(Anw);

[Anw, inds] = onlyGCC(Anw);
n = size(Anw,1);
comm1 = inds<=500;
A1 = Anw(comm1,comm1); A2 = Anw(~comm1,~comm1);
A12 = Anw(comm1,~comm1); A21 = Anw(~comm1,comm1);

deg = sum(Anw);
k = mean(deg); k2 = mean(deg.^2);
kappa = k2/k-1;

%% simulate

wOutVec = logspace(-1,1,50);
wInVec = logspace(0,2.8,50);
x = 0*wInVec'*wOutVec;
eta = x;

win_c1_mat = zeros(length(diff_conds),length(wOutVec));


for idx_conds = 1:length(diff_conds)
    conds = diff_conds{idx_conds};
    switch conds
        case 'free low'
            conditions.free_value = free_value_vec(1);
            NumExcited = NExcVec(1);
        case 'free high'
            conditions.free_value = free_value_vec(2);
            NumExcited = NExcVec(1);
        case 'ignited'
            conditions.free_value = free_value_vec(1);
            NumExcited = NExcVec(2);
    end
    
    for iwout = 1:length(wOutVec)
        wout = wOutVec (iwout)
        
        for iwin = 1:length(wInVec)
            win = wInVec (iwin);
            
            for ir = 1:reals
                s = randsample(find(comm1),1);
                conditions.specific_nodes = s;
                x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
                toplot = 0;
                                
                % find state 
                stateReals(ir) = find_num_comms_up_given_wout_win(win,wout);
                  
            end
            eta(iwin,iwout) = mean(stateReals);
        end
    end
end
toc


%% Figure rec vs w_intra w_inter
figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'layer','top','Xscale','log','yscale','log');
xlabel('\boldmath$\omega_{\rm intra}$','Interpreter','latex','FontSize',30)
ylabel('\boldmath$\omega_{\rm inter}$','Interpreter','latex','FontSize',30)

imagesc(wInVec,wOutVec,eta')
axis tight square
h = colorbar; h.Ticks = 0:2; h.LineWidth = 2;

xticks(10.^(-2:2:4)); yticks(10.^(-2:4))

urec_clr = [255 192 0]/255;
rec_clr = [0 100 50]/255;
rec2_clr = [0 51 102]/255;

alpha = linspace(0,1,50)';
map1 = (1-alpha)*urec_clr + (alpha)*rec_clr;
map2 = (1-alpha)*rec_clr + (alpha)*rec2_clr;
map = [map1; map2];

colormap(map)


%% save the figure in the folder 'output'
folder = '..\..\output\';
filename = '5g';
save_pdf_min_size([folder,filename])


%% examples

s = randsample(find(comm1),1);
conditions.specific_nodes = s;

%%

Delta = 5e3;
x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);

win = 20; wout = 5;

Xss = find_x_given_wout_win(win,wout,'all');

g = graph(Anw);
g.Nodes.activity = Xss;
g1 = subgraph(g,comm1); figure; p = plot(g1,'Layout','force','UseGravity','on'); x1 = p.XData/std(p.XData); y1 = p.YData/std(p.YData); close
g2 = subgraph(g,~comm1); figure; p = plot(g2,'Layout','force','UseGravity','on'); x2 = p.XData/std(p.XData); y2 = p.YData/std(p.YData); close
x2 = x2 + range(x1)*1.3;
x = [x1 x2]; y = [y1 y2];

figure; hold on
k = sum(A); [~, order] = sort(k); sizes = 2 + 4*order/n;
p = plot(g,'Layout','force','UseGravity','on','linewidth',1,'EdgeAlpha',0.25,'NodeCData',comm1,'MarkerSize',sizes);
p.XData = p.XData+0.3*comm1; p.YData = p.YData+0.3*comm1; p.NodeCData = double(comm1);
p.ZData = Xss;
axis square; colormap(flipud(jet))
view(40,20)

scatter3(p.XData(s),p.YData(s),p.ZData(s),130,'k','LineWidth',2,'MarkerFaceColor','y','MarkerFaceAlpha',0.5)
set(gca,'XTick',[],'YTick',[],'Visible','off','FontSize',20,'linewidth',2,'Box','on','ZScale','log');



%% save the figure in the folder 'output'
folder = '..\..\output\Figure5\';
filename = 'hij';
save_pdf_min_size([folder,filename])



%% functions
function x = find_x_given_wout_win(win,wout,all_or_eff)
global A A1 A2 A12 A21 x0 conditions M
A = [win*A1, wout*A12; wout*A21, win*A2];
x = SolveOdes(x0,M,all_or_eff, 0,conditions,0);
end

function num = find_num_comms_up_given_wout_win(win,wout)
global A A1 A2 A12 A21 x0 conditions M comm1
A = [win*A1, wout*A12; wout*A21, win*A2];
x = SolveOdes(x0,M,'all', 0,conditions,0);
x1 = mean(x(comm1)); x2 = mean(x(~comm1));
num = (x1>5) + (x2>5);
end


