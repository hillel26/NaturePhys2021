function movie_of_activity(t,x,x0)
global Anw

G = digraph(Anw);
s = find(x0);

G.Nodes.Activity = x(1,:)';
h = figure; hold off;
clim = [0,max(x(:))]; 
% clim = [0,11];
plot_Graph_Shells(G,s,h,clim)
zlim(clim);
set(gca,'nextplot','replacechildren');

m = length(t);
dm = ceil(m/10);
for i=1:dm:m
    G.Nodes.Activity = x(i,:)';
    h = figure;
    plot_Graph_Shells(G,s,h,clim);
end


