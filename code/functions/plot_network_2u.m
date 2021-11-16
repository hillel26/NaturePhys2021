function  [x,y] = plot_network_2u(g,coor)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    temp = figure;
    h_temp = plot(g,'layout','force');
    x = h_temp.XData; y = h_temp.YData;
    close(temp)
else
    x = g.Nodes.x;
    y = g.Nodes.y;
end

% edges
edges = table2array( g.Edges);
s = edges(:,1); t = edges(:,2);
m = length(s);
wmax = 2; % weight to compare to
widths =  g.Edges.Weight;
widths = 1 + (widths)/wmax*4;

map = [0 0.45 0.75; 0 0 0; 1 0 0];

figure; hold on
for i=1:m
    p1 = [x(s(i)) y(s(i))];
    p2 = [x(t(i)) y(t(i))];
    h_arcs(i) = arc(p1,p2,map(g.Edges.state(i),:),widths(i));
end

% nodes
n = g.numnodes;

NodesSizes = 50 + (g.degree-min(g.degree))/max(g.degree)*250;
NodesColors = g.Nodes.state;

scatter(x,y,NodesSizes,NodesColors,'o','filled','linewidth',2,'MarkerEdgeColor',0.2*[1 1 1])
% scatter(x(rmvnodes),y(rmvnodes),NodesSizes(rmvnodes),'k','o','filled','linewidth',2,'MarkerEdgeColor','k')
% scatter(x(prmvnodes),y(prmvnodes),NodesSizes(prmvnodes),[1 0 0],'o','filled','linewidth',2,'MarkerEdgeColor',0.6*[1 0 0])

colormap(map)
axis equal
set(gca,'visible','off')
set(gcf,'Color','w')


