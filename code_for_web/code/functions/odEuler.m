function [t,x] = odEuler(eq,time,x0)
dt=1e-6;

tt = time(1);
xt = x0;

t = zeros(50,1);
x = zeros(length(t),length(x0));

it = 1;
t(it) = tt;
x(it,:) = xt;

Di = floor(range(time)/length(t)/dt);
i=0;
while tt<time(end)
    i=i+1;
    xt = xt + dt*eq(t,xt);
    tt = tt+dt;
    if ~mod(i,Di)
        it=it+1;
        t(it) = tt;
        x(it,:) = xt;
    end
end

end
