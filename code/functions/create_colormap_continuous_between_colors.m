function map = create_colormap_continuous_between_colors(c1,c2,c3)
a = linspace(0,1,51)';
map = a*c2+(1-a)*c1;
if nargin>2
    
    a1 = 1-a;
    a3 = a;
    a21 = 2*a(1:ceil(end/2)).^2;
    a22 = flipud(a21);
    a22(1) = [];
    a2 = [ a21; a22];
    
    map = a1*c1 + a2*c2 + a3*c3;
    map = map./max(map,[],2);
    
%     map1 = a.^3*c2+(1-a.^3)*c1;
%     a(1)=[];
%     map2 = flipud( a.^2*c2+(1-a.^2)*c3 );
%     map = [map1;map2];
end

