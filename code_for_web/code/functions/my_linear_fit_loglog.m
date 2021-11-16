function [p] = my_linear_fit_loglog(x,y,side,width,color,name_of_exponent)

% default settings
if nargin<3; side = 'right'; end
if nargin<4; width = 1/5; end
if nargin<5; color = 'k'; end
if nargin>5
    name = [name_of_exponent,'='];
else
    name = '';
end

p = polyfit(log10(x),log10(y),1);
fit = @(x) 10.^polyval(p,log10(x));
hold on
plot(x,fit(x),'color',color,'Linewidth',1)
points = sqrt(x(1)*x(end))*(x(end)/x(1)).^([-1,1]*width/2);
points = sort(points);

switch side
    case 'right'
        plot(repelem(points,[1,2]), ...
            fit( repelem(points,[2,1]) )...
            ,'color',color,'linewidth',1);
        text(points(2), fit(geomean(points)),...
            [name,num2str(p(1),2)],'fontsize',15,'color',color,...
            'horizontalalignment','left','VerticalAlignment','middle')
    case 'left'
        plot(repelem(points,[2,1]), ...
            fit( repelem(points,[1,2]) )...
            ,'color',color,'linewidth',1);
        text(points(1),fit(geomean(points)),...
            [name,num2str(p(1),2)], 'fontsize',15,'color',color,...
            'horizontalalignment','right','VerticalAlignment','middle')
end

set(gca,'Xscale','log','Yscale','log')
