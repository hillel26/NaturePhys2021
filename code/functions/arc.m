function h_arc = arc(p1,p2,color,linewidth)
rng(0)
beta = pi/4;
s = p2-p1;
n = [s(2),-s(1)]/norm(s) * (-1)^randi(2);
o = p1+s/2 + n*norm(s)/2*cot(beta/2);
r = norm(p1-o);
th1 = acos((p1(1)-o(1))/r);
if p1(2)<o(2); th1 = -th1; end
th2 = acos((p2(1)-o(1))/r);
if p2(2)<o(2); th2 = -th2; end
th12 = sort([th1,th2]);
if diff(th12) > pi; th12(1) = th12(1) + 2*pi; end
th = linspace(th12(1),th12(2));
x = o(1) + r*cos(th);
y = o(2) + r*sin(th);
h_arc = plot(x,y,'color',color,'linewidth',linewidth);
end