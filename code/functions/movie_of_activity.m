function movie_of_activity(t,x,x0)
global Anw

G = graph(Anw);
s = find(x0);

G.Nodes.Activity = x(1,:)';
h = figure; hold off;
clim = [0,max(x(:))]; 
% clim = [0,11];
plot_Graph_Shells(G,s,h,clim)
zlim(clim);
set(gca,'nextplot','replacechildren');

randnum = num2str(randi(1000));
name = [ 'illustrate',randnum,'.avi'];
v = VideoWriter(name);
v.FrameRate = 10;
v.Quality = 100; 
open(v);

m = length(t);
dm = ceil(m/400);
for i=1:dm:m
    G.Nodes.Activity = x(i,:)';
    plot_Graph_Shells(G,s,h,clim);
    axis off
%     text(1/2,1,['t=',num2str(t(i))],'Units','normalized','HorizontalAlignment','center','FontSize',15)
    frame = getframe(gcf);
    writeVideo(v,frame);
%     if ~ mod(i-1,5)
%         saveas(gcf,[randnum,'frame',num2str(i)],'fig')
%     end
end

close(v);
