%% Comparison between tree and ignite
%% MicroBiomeDynamics

%% introductions
clear; clc;  close all; format compact;

addpath('..\functions')

active_color = [0 210 200]/255; %[0 170 170]/255;
bi_color = 0.3*[1 1 1];
inactive_color = 'r'; %[155 0 0]/255;

global A Anw x0  NameOfModel Factor_Xeff F B

NameOfModel = 'Eco'; %
F = 5; % parameter of Eco
B = 3; % parameter of Eco
M = KindOfDynamics( NameOfModel );

load('..\..\data\MicrobiomeNetworks')

ChooseExcNodes = 'specific'; % % options - {'specific', 'low degree', 'high degree', 'local', 'rand'}
conditions.type =  'BC'; % 'pills'; %
conditions.TimeProBio = 5;
all_or_eff = 'average'; % what xss we need: 'all'/ 'eff' /'eff_F'
release = 0; % to release or not after holding

%% build the matrix
n = size(competition,1);
Ap = 30; Ac = 1; % positive weight (wp) and competiton weight (wc)
P = Ap*complementarity; % matrix of positive dependency interactions
C = Ac*competition; % matrix of competitive interactions

A = sparse(P-C);

% plot A
Clrs = sign(A);
figure; hold on
imagesc(Clrs); axis tight square off
colormap([0.8 0.3 0.1; 1 1 1; 0 0.45 0.75])

% save the figure in the folder 'output'
folder = '..\..\output\Figure6\';
filename = 'c';
save_pdf_min_size([folder,filename])

%% find the upper state
conditions.free_value = 10;
conditions.Delta = 0;
toplot = 0;
Factor_Xeff = sum(A,1)/sum(sum(A));
[Xss] = SolveOdesProbiotic(zeros(n,1), M, 'all', toplot ,conditions,release);
% figure; hold on; histogram(Xss,0:35,'FaceColor',active_color); histogram(Xss,0:1,'FaceColor',inactive_color)
% set(gca,'Fontsize',15,'box','on','LineWidth',2); xlabel('$x$','Interpreter','latex','FontSize',20); ylabel('Frequency'); 

%% remove nodes with low activity
active_nodes = find(Xss>2);
A = A(active_nodes,active_nodes);
Anw = logical(A);
n = size(A,1);

% plot A
Clrs = sign(A);
figure; hold on
imagesc(Clrs); axis tight square off
colormap([0.8 0.3 0.1; 1 1 1; 0 0.45 0.75])

% save the figure in the folder 'output'
filename = 'o_lu';
save_pdf_min_size([folder,filename])

close_up = 50;
figure; hold on
imagesc(Clrs); axis tight square off
colormap([0.8 0.3 0.1; 1 1 1; 0 0.45 0.75])
axis(-0.5 + [1 close_up 1 close_up])

% save the figure in the folder 'output'
filename = 'o_ld';
save_pdf_min_size([folder,filename])

%% find now the upper state
conditions.free_value = 10;
conditions.Delta = 0;
Factor_Xeff = sum(A,1)/sum(sum(A));
[Xss] = SolveOdesProbiotic(zeros(n,1), M, 'all', toplot ,conditions,release);
% figure; hold on; histogram(Xss,0:35,'FaceColor',active_color); histogram(Xss,0:1,'FaceColor',inactive_color)
% set(gca,'Fontsize',15,'box','on','LineWidth',2); xlabel('$x$','Interpreter','latex','FontSize',20); ylabel('Frequency'); 

% plot the state
G = digraph(A');
h = figure;
clim = [0, 40];
G.Nodes.Activity = Xss;
plot_Graph_Shells(G,[],h,clim);
p = h.Children.Children(1);
pos_edges = G.Edges.Weight'>0; neg_edges = G.Edges.Weight'<0;
p.EdgeAlpha= 0.02;  
highlight( p, 'Edges',neg_edges,'EdgeColor',[0.8 0.3 0.1],'LineWidth',0.1);

% save the figure in the folder 'output'
filename = 'e';
save_pdf_min_size([folder,filename])


%% find now the lower state
conditions.free_value = 0.1;
conditions.Delta = 0;
Factor_Xeff = sum(A,1)/sum(sum(A));
[Xss] = SolveOdesProbiotic(zeros(n,1), M, 'all', 0 ,conditions,release);
% figure; hold on; histogram(Xss,0:35,'FaceColor',active_color); histogram(Xss,0:1,'FaceColor',inactive_color)
% set(gca,'Fontsize',15,'box','on','LineWidth',2); xlabel('$x$','Interpreter','latex','FontSize',20); ylabel('Frequency'); 

% plot the state
G = digraph(A');
h = figure;
clim = [0, 40];
G.Nodes.Activity = Xss;
plot_Graph_Shells(G,[],h,clim);
p = h.Children.Children(1);
pos_edges = G.Edges.Weight'>0; neg_edges = G.Edges.Weight'<0;
p.EdgeAlpha= 0.02;  
highlight( p, 'Edges',neg_edges,'EdgeColor',[0.8 0.3 0.1],'LineWidth',0.1);

% save the figure in the folder 'output'
filename = 'f';
save_pdf_min_size([folder,filename])


%% simulate reviving by each node

conditions.free_value = 0.1;
NumExcited = 1;

Delta = 10;
conditions.Delta = 10;
toplot = 0;

Factor_Xeff = sum(A,1)/sum(sum(A));

state = zeros(n,1);

for i=1:n
    i
    conditions.specific_nodes = i; % find(active_nodes==189); % 140 129 673 157 %ignite_nodes(30);
    
    % initial/fixed condition
    x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
    
    state(i) = SolveOdesProbiotic(x0, M, all_or_eff, toplot ,conditions,release);
end


%% find reaches and rank all species
G = digraph(A');

w0 = 1.122;
g = G.rmedge(find(G.Edges.Weight<w0)); % w_th for a tree spreading the reiginition
Reaches = arrayfun(@(s) length( g.bfsearch(s) ) ,(1:n)');
d = g.distances('Method','unweighted'); d(isinf(d))= nan; d = max(d',[],'omitnan')';

% ranking species

[~, ranked_nodes] = sort(Reaches,'descend');

% compare to simulations
figure; hold on
ylim([-0.1 1.1]); 
patch([0 27 27 1],repelem(ylim,2),0.75*[1 1 1],'EdgeColor','none')
scatter(1:n,state(ranked_nodes)>5,10,'filled')
xlabel('Rank','FontSize',20); ylabel('\boldmath$\eta$','Interpreter','latex','FontSize',30)
set(gca,'Fontsize',20,'linewidth',2,'box','on','XScale','lin','yscale','lin','layer','top')
yticks([0 1])
xlim([-20 n+20]);

% save the figure in the folder 'output'
folderSI = '..\..\output\FiguresSI\';
filename = '7';
save_pdf_min_size([folderSI,filename])

% save in table
Rank = (1:n)';
IgniterSpeciesIndex = active_nodes(ranked_nodes(Rank));
IgniterSpecies = species(IgniterSpeciesIndex);
TreeSize = Reaches(ranked_nodes(Rank));
% Shells = d(ranked_nodes(Rank));
State = 1/10*round(10*state(ranked_nodes(Rank)));

T1 = table( Rank, IgniterSpecies, TreeSize, State);
writetable(T1,[folder,'i.xlsx']);
writetable(T1,[folderSI,'dataset1.xlsx']);



%% plot reviving and reach for speciel species
conditions.free_value = 0.1;
NumExcited = 1;

pnls = 'jkl';
ip = 0;

for igniter = [116, 642 582] %  116 - Bifidobacterium bifidum ; 642 - Pseudomonas putida ; 582 - Parabacteroides distasonis 
conditions.specific_nodes = find(active_nodes==igniter); 
species(igniter)


Factor_Xeff = sum(A,1)/sum(sum(A));

x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);

[Xss] = SolveOdes(x0, M, 'all', toplot ,conditions,release);

% plot reach
s = conditions.specific_nodes;
H = g.subgraph(g.bfsearch(s));
H.Nodes.Activity = Xss(g.bfsearch(s));
clim = [0, 40];
h = figure;
plot_Graph_Shells(H,1,h,clim)
p1 = h.Children.Children(2);
view(2)
axis equal off

ip = ip+1;

% save the figure in the folder 'output'
filename = [pnls(ip),'1'];
save_pdf_min_size([folder,filename])

% plot state
h = figure;
G.Nodes.Activity = Xss;
plot_Graph_Shells(G,s,h,clim);
p = h.Children.Children(2);
pos_edges = G.Edges.Weight'>0; neg_edges = G.Edges.Weight'<0;
p.EdgeAlpha= 0.02;  
highlight( p, 'Edges',neg_edges,'EdgeColor',[0.8 0.3 0.1],'LineWidth',0.1);

% save the figure in the folder 'output'
filename = [pnls(ip),'2'];
save_pdf_min_size([folder,filename])

end

%% rank metabolites impact on the competition
m = size(norm_import,2);

idxs = 1:m;
Q = norm_import*norm_import';
QT =sum(Q(:));
impact = zeros(m,1);
for i=1:m
    Qi = norm_import(:,idxs~=i)*norm_import(:,idxs~=i)';
    QiT = sum(Qi(:));
    impact(i) = 1-QiT/QT;
end 

[SortImpact, SortMetabolites] = sort(impact,'descend');
Rank = (1:m)';
TopMetabolites = metabolites(SortMetabolites(Rank));
TopImpact = SortImpact(Rank);
T2 = table( Rank, TopMetabolites, TopImpact);
writetable(T2,[folder,'n.xlsx'])
writetable(T2,[folderSI,'dataset2.xlsx']);

%% remove top metabolites from competition
rmv_metab = SortMetabolites(1:3);
norm_import(:,rmv_metab) = 0;
competition = norm_import*norm_import';
competition(1:n+1:end) = 0; % diagonal is zero
C = Ac*competition; 

A = sparse(P-C);
A = A(active_nodes,active_nodes);

% plot A
Clrs = sign(A);
figure; hold on
imagesc(Clrs); axis tight square off
colormap([0.8 0.3 0.1; 1 1 1; 0 0.45 0.75])

% save the figure in the folder 'output'
filename = 'o_ru';
save_pdf_min_size([folder,filename])

close_up = 50;
figure; hold on
imagesc(Clrs); axis tight square off
colormap([0.8 0.3 0.1; 1 1 1; 0 0.45 0.75])
axis(-0.5 + [1 close_up 1 close_up])

% save the figure in the folder 'output'
filename = 'o_rd';
save_pdf_min_size([folder,filename])


%% plot reviving and reach for speciel nodes after removal of top competitions
conditions.free_value = 0.1;
NumExcited = 1;

igniter = 642; % Pseudomonas putida
conditions.specific_nodes = find(active_nodes==igniter); 
species(igniter)

Factor_Xeff = sum(A,1)/sum(sum(A));

x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);

[Xss] = SolveOdes(x0, M, 'all', toplot ,conditions,release);

% plot state
s = conditions.specific_nodes;
h = figure;
G.Nodes.Activity = Xss;
plot_Graph_Shells(G,s,h,clim);
p = h.Children.Children(2);
pos_edges = G.Edges.Weight'>0; neg_edges = G.Edges.Weight'<0;
p.EdgeAlpha= 0.02;  
highlight( p, 'Edges',neg_edges,'EdgeColor',[0.8 0.3 0.1],'LineWidth',0.1);

% save the figure in the folder 'output'
filename = 'q';
save_pdf_min_size([folder,filename])




