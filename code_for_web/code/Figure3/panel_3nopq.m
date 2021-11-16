%% introductions
clear; clc; format compact; % close all

addpath('..\functions');

global Anw x0 NameOfModel
tic

%% Dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
% a=1; b=2; % parameters of MM
M = KindOfDynamics( NameOfModel );
range_sol = [0,10]; % range of solutions for high and low states in x

%% Reviving
conditions.type = 'BC';
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand'}
Delta = 1; % 10; % force of reigniting - holding value

all_or_eff = 'eff'; % what xss we need: 'all'/ 'eff' /'eff_F'
free_value = 0; % the initial value of the free nodes
NumExcited = 1; %round(N*[0,0.1]); % round(N*logspace(-4,-1,30)); %
release = 1; % to release or not after holding
reals = 50; % number of realizations
th = 1/2; % threshold determining seccuess of reviving (if Xeff>th)

conditions.free_value = free_value;

%% The network
NetStruct =  'PPI_Yeast'; % 'PPI_Human'; %  'ER' ; %  'SF'; %  
N = 1e4; k = 10; gamma = 2.7; k0 = 2.6;
Anw = BuildNetwork(N, NetStruct,k,'gcc'); % Adjacency matrix, not weighted
n = size(Anw,1);
% k = mean(sum(Anw));
% kp2 = mean(sum(Anw).^2);
% kappa = kp2/k-1;

kin = sum(Anw,2);
kout = sum(Anw,1)';
kappa = mean( (Anw'*(kin-1))./kout ); % according to the definition in the paper

%% simulations
wVec = logspace(-1.5,1,50);
success_rate = 0*wVec;


for iw=1:length(wVec)
    w = wVec(iw);
    disp(iw/length(wVec))
    
    successes = 0;
    for ir = 1:reals
        
        % initial/fixed condition
        x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,Delta);
        toplot=0;
                
        Xss = find_x_via_omega(w,M,all_or_eff,toplot,conditions,release);
        successes = successes + (Xss>th);
    end
    success_rate(iw) = successes/reals;
end

toc

%% Theory for recoverability
options = optimset('TolX',1e-3);
wc = fzero(@(w) is_rcoverable_theory(w,kappa,Delta,range_sol) - th , [0.1, 2],options);

%% figure

figure; hold on
scatter(wVec,success_rate,80,'o','LineWidth',2,'MarkerEdgeColor',[0 51 102]/255) % [0 0.7 1]
set(gca,'FontSize',25,'linewidth',2,'box','on','XScale','log','YScale','lin')
xlabel('\boldmath$\omega$','Interpreter','latex','FontSize',35);
ylabel('\boldmath$\eta$','Interpreter','latex','FontSize',35);
axis tight
xlim(10.^[-1.5, 1]); xticks(10.^(-2:2)); ylim([-0.05,1.05])
plot(wc*[1 1],ylim,'--','Color',[0 0 0],'LineWidth',2)

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'nopq';
save_pdf_min_size([folder,filename])
