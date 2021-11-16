clear; clc; close all
global A Anw x0 NameOfModel

addpath('..\functions');

NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
M = KindOfDynamics( NameOfModel );

conditions.type = 'BC';
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
release = 1; % to release or not after holding6

conditions.free_value = 10; % the initial value of the free nodes
NumExcited = 0; %round(N*[0,0.1]); % round(N*logspace(-4,-1,30)); %
Delta = 5;
range = [0,20];

%% get the network

NetStruct =  'PPI_Yeast'; % 'PPI_Human'; %  'ER' ; %  'SF'; %
N = 0; k = 0; % not relevant but needs some input
Anw = BuildNetwork(N, NetStruct, k, 'gcc'); % Adjacency matrix, not weighted

A_original = Anw;

g_or = graph(A_original);
edges = table2array( g_or.Edges);
s = edges(:,1); t = edges(:,2);

n = g_or.numnodes;
m = g_or.numedges;

kappa1 = mean(sum(Anw).^2)/mean(sum(Anw))-1

%% plot the complete network and activity
figure
plot(g_or,'Layout','force')
axis off

x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
g_or.Nodes.Activity = find_x_via_omega(1,M,'all',0,conditions,release);
clim = [0 max(g_or.Nodes.Activity)+5];
plot_Graph_Shells(g_or,[],figure,clim )


%% damage towards the inactive phase

qlinks = 0.3;
rmvlinks = randsample(m,round(qlinks*m));
for i=1:length(rmvlinks)
    Anw(s(rmvlinks(i)),t(rmvlinks(i))) = 0;
    Anw(t(rmvlinks(i)),s(rmvlinks(i))) = 0;
end

qnodes = 0.3;
rmvnodes = randsample(n, round(qnodes*n)); % removed nodes
Anw(rmvnodes,:) = 0; % nodes loss
Anw(:,rmvnodes) = 0; % nodes loss

% Anw = onlyGCC(Anw);

k = mean(sum(Anw))
kappa0 = mean(sum(Anw).^2)/mean(sum(Anw))-1
w0 = 0.2

A = 0.2*Anw;


%% constraints

qplinks = 0.5; % #  permanent removed links
prmvlinks = randsample(rmvlinks, round(qplinks*length(rmvlinks))); % permanent removed links

qpnodes = 0.5; % #  permanent removed nodes
prmvnodes = randsample(rmvnodes, round(qpnodes*length(rmvnodes))); % permanent removed nodes

g = graph(A);
g_or.Edges.Weight(:) = w0;


%% plot the network
g_or.Nodes.state(:) = 1;
g_or.Nodes.state(rmvnodes) = 2;
g_or.Nodes.state(prmvnodes) = 3;
g.Nodes.state(:) = 1;

g_or.Edges.state(:) = 1;
g_or.Edges.state(rmvlinks) = 2;
g_or.Edges.state(prmvlinks) = 3;
g.Edges.state(:) = 1;

figure
p = plot(g_or,'Layout','force');

p.EdgeCData = g_or.Edges.state;
p.NodeCData = g_or.Nodes.state;

map = [0 0.45 0.74; 0 0 0; 1 0 0];
colormap(map)

set(gca,'Visible','off')


%% calculate and plot the state

g_or.Nodes.Activity = find_x_via_omega(w0,M,'all',0,conditions,release);
clim = [0 max(g_or.Nodes.Activity)+5];
plot_Graph_Shells(g_or,[],figure,clim )



%% extract small sub graph and plot

[bins, binsize] = conncomp(g);
[~,gccnum] = max(binsize);
gcc = find(bins==gccnum);

root = randsample(gcc,1);
v2 = bfsearch(g,root);
v2 = v2(1:min(10,end));
v1 = bfsearch(g_or,root);
for l=1:length(v1)
    if all(ismember(v2,v1(1:l)))
        break
    end
end
v = v1(1:l+15);

gsub1 = subgraph(g_or,v);

[xv, yv] = plot_network_2u(gsub1);
g.Nodes.x(v) = xv; g.Nodes.y(v) = yv;

gsub2 = subgraph(g, v);
[bins, binsize] = conncomp(gsub2);
[~,gccnum] = max(binsize);
gcc2 = find(bins==gccnum);
gsub2 = subgraph(gsub2, gcc2);

plot_network_2u(gsub2,'use_coordinates');
colormap([0 0.45 0.75])

% g_failed = g;
A_failed = A;

%% gain weight nodes and links

fig = open('BackgroundForPanel4f.fig'); hold on; ax = fig.Children;

ip = 0;
for p = [0.4 0.6 0.9]
    ip=ip+1
    % g = g_failed; % in order to return to the same network
    A = A_failed;
    rmvnodes = find(g_or.Nodes.state~=1);
    rmvlinks = find(g_or.Edges.state~=1);
    w = [];
    kappa = [];
    
    for idx_build=1:1e4
        if rand<p
            for j=1:5
                increaseWeight = randsample(setdiff(1:n,rmvlinks),1);
                if all(not( ismember([s(increaseWeight),t(increaseWeight)], rmvnodes)))
                    A(s(increaseWeight),t(increaseWeight)) = 3*w0*rand + A(s(increaseWeight),t(increaseWeight));
                    A(t(increaseWeight),s(increaseWeight)) = A(s(increaseWeight),t(increaseWeight));
                end
            end
        elseif rand<0.25 && not(isempty(setdiff(rmvnodes,prmvnodes)))
            returnnode = randsample(setdiff(rmvnodes,prmvnodes),1);
            A(returnnode, setdiff( g_or.neighbors(returnnode), rmvnodes)) = w0;
            A(setdiff( g_or.neighbors(returnnode), rmvnodes), returnnode) = w0;
            %         g_or.Nodes.state(returnnode) = 1;
            rmvnodes = setdiff(rmvnodes,returnnode);
        elseif not(isempty( setdiff(rmvlinks,prmvlinks)))
            returnlink = randsample(setdiff(rmvlinks,prmvlinks),1);
            A(s(returnlink), t(returnlink)) = w0;
            A(t(returnlink), s(returnlink)) = w0;
            %     g_or.Edges.state(returnlink) = 1;
            rmvlinks = setdiff(rmvlinks,returnlink);
        end
        
        kappa(idx_build) = mean(sum(logical(A)).^2)/mean(sum(logical(A))) - 1;
        w(idx_build) = mean(A(A(:)~=0));
        
        if ip == 2 && ~mod(idx_build,400)
            g = graph (A);
            g.Nodes.x(v) = xv; g.Nodes.y(v) = yv;
            
            g.Nodes.state(:) = 1; g.Edges.state(:) = 1;
            gsub2 = subgraph(g, v);
            [bins, binsize] = conncomp(gsub2);
            [~,gccnum] = max(binsize);
            gcc2 = find(bins==gccnum);
            gsub2 = subgraph(gsub2, gcc2);
            
            plot_network_2u(gsub2,'use_coordinates');
            colormap([0 0.45 0.75])
        end
        
        % last figure and home
        if w(idx_build)> 0.95 - 0.02*kappa(idx_build)
            g = graph (A);
            g.Nodes.x(v) = xv; g.Nodes.y(v) = yv;
            
            g.Nodes.state(:) = 1; g.Edges.state(:) = 1;
            gsub2 = subgraph(g, v);
            [bins, binsize] = conncomp(gsub2);
            [~,gccnum] = max(binsize);
            gcc2 = find(bins==gccnum);
            gsub2 = subgraph(gsub2, gcc2);
            
            plot_network_2u(gsub2,'use_coordinates');
            colormap([0 0.45 0.75])
            break
        end
    end
    
    if ip==2
        figure(100); hold on
        set(gca,'box','on','FontSize',20,'LineWidth',2)
        yyaxis('left')
        ylabel('\boldmath$\kappa$','Interpreter','latex','FontSize',30)
        plot(kappa,'LineWidth',2)
        yyaxis('right')
        ylabel('\boldmath$w$','Interpreter','latex','FontSize',30)
        plot(w,'LineWidth',2)
        axis tight
    end
    
    plot(ax,kappa,w,'LineWidth',2)
    
    
    % calculate and plot the state
    NumExcited = 1;
    Delta = 10;
    conditions.free_value = 0;
    release = 0;
    A = onlyGCC(A);
    Anw = logical(A);
    n = size(A,1);
    x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
    
    SolveOdes(x0,M,'all','pics',conditions,release);
    clim = [0 max(g_or.Nodes.Activity)+5];
    % plot_Graph_Shells(g_or,[],figure,clim )
    
end


%% save only figure 4e and 4f in the folder 'output'
folder = '..\..\output\Figure4\';
filename = 'e';
figure(100)
save_pdf_min_size([folder,filename])

filename = 'f';
figure(fig.Number)
save_pdf_min_size([folder,filename])



%% diagram

% kappa = logspace(0.5,1.5,15);
% options = optimset('TolX',1e-3);
% wc_theory = arrayfun(@(kappa) fzero(@(w) is_rcoverable_theory(w,kappa,holding_value,range) - 1/2 , [2/(kappa+1),1e1],options) , kappa);
%
% % inactive phase
% wc_inactive_theory = 2./(kappa+1);
%
% figure; hold on
% set(gca,'FontSize',20,'box','on','LineWidth',1.5,'XScale','log','YScale','log','layer','top')
% xlabel('\boldmath$\kappa$','Interpreter','latex','FontSize',33)
% ylabel('\boldmath$w$','Interpreter','latex','FontSize',33)
% axis square
%
% ymin = 10^-1.5; ymax = 4;
% xlim(kappa([1,end]))
% ylim([ymin,ymax])
%
% % patches
% clrs = [0 51 102; 255 192 0; 155 0 0]/255;
%
% patch([xlim,fliplr(xlim)],repelem(ylim,2),clrs(2,:))
% patch([kappa,fliplr(kappa)],[ymax+0*kappa,fliplr(wc_theory)],clrs(1,:))
% patch([kappa,fliplr(kappa)],[ymin+0*kappa,fliplr(wc_inactive_theory)],clrs(3,:))
%
%
% % line and symbols
% % loglog(kappa,wc_theory,'w','LineWidth',2)
% % loglog(kappa,wc_inactive_theory,'w','LineWidth',2)
%
% xticks(10:10:30)

