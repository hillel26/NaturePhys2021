function [ Xss ,tau1, tau2] = SolveOdes(x0, M, all_or_eff, toplot ,conditions,release,times)
%SOLVEODES - solves set of ode's
%   get (A M0 M1 M2 omega x0) and give Xss
%   all_or_eff = 'all' - vector of xi eff / 'eff' - scalar of Xeff /...
%   'eff_F' - scalar of Xeff_F only on the free set F.
%   toplot = 1 to plot xi_vs_t / 0 not
%   conditions_type is 'IC' or 'BC' initial or boundary conditions
%   'times' is if you want to calculate relaxation times and how - {'half-life','steady-state'}

global A Anw Factor_Xeff NameOfModel

M0 = M{1};
M1 = M{2};
M2 = M{3};

switch conditions.type
    case 'IC'
        Equations = @(t,x) (M0(x) + M1(x).*(A*M2(x)));
        x_free = not(logical(x0));
        x0(x_free) = conditions.free_value;
    case 'BC'
        x_free = not(logical(x0));
        Equations = @(t,x) (M0(x) + M1(x).*(A*M2(x))).*x_free;
        if strcmp(NameOfModel,'Glauber')
            Equations = @(t,x) (M0(x) + M1(x).*(M2(A*x))).*x_free;
        elseif strcmp(NameOfModel,'Ising-Sch')
            Equations = @(t,x) (M0(x) + A*M2(x)-x.*sum(A,2)).*x_free;
        end
        x0(x_free) = conditions.free_value;
end

% solve the ode!!
T = 0;
dT = 5;
Tmax = 20; % !!!!!!!!!!!!!!!
t = [];
x = [];
delta = 1; % in order to run first iteration
while   T<Tmax && delta > 1e-3
%     options = odeset('reltol',1e-8,'abstol',1e-10);
    [t1,x1] = ode45(Equations,linspace(T,T+dT,100),x0); % 'options
    t = [t;t1];
    x = [x;x1];
    [delta, ~] = max(abs(( mean(x(end-5+1:end,:)) - mean(x(end-10+1:end-5,:)) )/( t(end-2)-t(end-7) )));
    T = T+dT;
    x0 = x(end,:)';
end
idx_ss1 = length(t); % using for the relaxation time

% release the fixed nodes

% x0(~x_free)=0; % !!!!!!!!!!!!!!!!!!!!!!!!
Tmax = Tmax + T-dT;
if strcmp(conditions.type,'BC') && release
    Equations = @(t,x) (M0(x) + M1(x).*(A*M2(x)));
    delta = 1; % in order to run first iteration
    while   T<=Tmax && delta > 1e-2
        options = odeset('reltol',1e-5,'abstol',1e-8);
        [t1,x1] = ode45(Equations,[T T+dT],x0,options);
        t = [t;t1];
        x = [x;x1];
        [delta, ~] = max(abs(( mean(x(end-5+1:end,:)) - mean(x(end-10+1:end-5,:)) )/( t(end-2)-t(end-7) )));
        T = T+dT;
        x0 = x(end,:)';
    end
end

% output
switch all_or_eff
    case 'all'
        Xss = x0;
    case 'eff'
        Xss = Factor_Xeff * x0; % Xeff
    case 'eff_F'
        Factor_Xeff_F = sum(A(x_free,x_free),1)/sum(sum(A(x_free,x_free)));
        Xss = Factor_Xeff_F * x0(x_free); % Xeff
    case 'average'
        Xss = mean(x0);
    otherwise
        disp('error , first arg should be ''all'' or ''eff''')
end

% plot
if strcmp(toplot,'all') || strcmp(toplot,'time') 
    figure
    subplot(2,1,1)
    plot(t,x)
    xlabel('t')
    ylabel('x')
    subplot(2,1,2)
    plot(t,x*Factor_Xeff')
    xlabel('t')
    ylabel('x_{eff}')
    w = sum(A(:))/nnz(A); k = sum(Anw); kappa = mean(k.^2)/mean(k); beta = w*kappa;
    title(['w=',num2str(w),' \kappa=',num2str(kappa),' \beta=',num2str(beta),])
end
if strcmp(toplot,'all') || strcmp(toplot,'movie') 
    movie_of_activity(t,x,~x_free)
end
if strcmp(toplot,'all') || strcmp(toplot,'pics') 
    pics_of_activity(t,x,~x_free)
end

% calculate the relaxation times
if nargin>6
    switch times
        case 'half-life'
            alpha = 0.9;
            beta = 0.9;
            xeff_vs_t = x*Factor_Xeff';
            tau1 = t(find( abs(xeff_vs_t-xeff_vs_t(1)) > alpha * abs(xeff_vs_t(idx_ss1)-xeff_vs_t(1)) , 1 ) );
            tau2 = t(idx_ss1 + ...
                find( abs(xeff_vs_t((idx_ss1+1):end)-xeff_vs_t(idx_ss1+1)) > ...
                beta * abs(xeff_vs_t(end)-xeff_vs_t(idx_ss1+1)) , 1 ) )...
                - t(idx_ss1);
        case 'steady-state'
            tau1 = t(idx_ss1);
            tau2 = t(end)-tau1;
    end
end
end

