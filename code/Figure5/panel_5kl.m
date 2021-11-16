%% w_inter vs w_intra diagram

clear; format compact;
close all;

global A Anw x0 A1 A2 A12 A21 M comm1 mu

addpath('..\functions');

load('..\..\data\Brain')

tic

%% dynamics
NameOfModel = 'Neural'; % 'Eco'; % 'MM'; %  'Simple'; % 'Glauber'; % 'Voter'; % 'SIS';% 'MAK';%  'PD';%
mu = 10;

M = KindOfDynamics( NameOfModel );
range_sols = [0,50]; % range of high and low state for the search of the theory solution

%% build the network 
% remove weak links

rng(0)

Anw = (A>0.03);
Anw = sparse(Anw);

[Anw, inds] = onlyGCC(Anw);
n = size(Anw,1);
comm1 = inds<=500;
A1 = Anw(comm1,comm1); A2 = Anw(~comm1,~comm1);
A12 = Anw(comm1,~comm1); A21 = Anw(~comm1,comm1);

degs = sum(Anw);
k = mean(degs); k2 = mean(degs.^2);
kappa = k2/k-1;

% set different weights
wout = 5;
win = 20;

A = [win*A1, wout*A12; wout*A21, win*A2];

%% simulate 

M0 = M{1}; M1 = M{2}; M2 = M{3};
Equations = @(t,x) (M0(x) + M1(x).*(A*M2(x)));
x0 = 30*rand(n,1);

T = 0;
dT = 1;
Tmax = 100; 
t = []; x1 = []; x2 = [];

while  T<Tmax
    [tnew,x] = ode45(Equations,linspace(T,T+dT,100),x0); % 'options
    t = [t;tnew];
    x1new = mean(x(:,comm1),2); x2new = mean(x(:,~comm1),2);
    x1 = [x1;x1new];  x2 = [x2;x2new]; 

    T = T+dT;
    x0 = x(end,:)';
    
    % noise
    x0(comm1) = x0(comm1) * rand;
    x0(~comm1) = x0(~comm1) * rand;
end                            
                              
toc

%% Figure rec vs w_intra w_inter
figure; hold on
set(gca,'FontSize',20,'box','on','LineWidth',2,'layer','top','Xscale','lin','yscale','log');
ylabel('\boldmath$\bar{\rm x}$','Interpreter','latex','FontSize',30)
xlabel('\boldmath$t$','Interpreter','latex','FontSize',30)
set(gcf,'Position',[200 200 650 300])

colororder = [0 0 0.51; 0.5 0 0];

plot(t,x1,'LineWidth',2,'Color',[0 0 0.51])
plot(t,x2,'LineWidth',2,'Color',[0.5 0 0])

ax = gca;
ylim([1e-1,1e2]); % ylim([1e-5,1e2])
ax.YTick = 10.^(-4:2:3);


%% save the figure in the folder 'output'
folder = '..\..\output\Figure5\';
filename = 'kl';
save_pdf_min_size([folder,filename])

