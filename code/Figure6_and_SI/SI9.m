clear; clc; % close all
format compact;

addpath('..\functions');

global A Anw x0 Factor_Xeff List lambda NameOfModel a b
tic

% model of dynamics
NameOfModel = 'MM'; %  'Neural'; % 'Eco'; %   'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%
% a=1; b=2; % parameters of MM
M = KindOfDynamics( NameOfModel );
conditions.type =  'pills'; %  'IC'; %    'BC'; % 
ChooseExcNodes = 'rand'; % options - {'low degree', 'high degree', 'local', 'rand','specific'}
holding_value = 2;
conditions.Delta = 2;
TimeProBioVec = 2:4;% 10; % 


all_or_eff = 'eff_F'; % what xss we need: 'all'/ 'eff' /'eff_F'
release = 1; % to release or not after holding
conditions.free_value = 0;
NumExcited = 1;

N=1e4;

NetStruct = 'ER';
k = 5; %round(logspace(log10(10),log10(20),4));
w = 0.8; % weight

%% build the network
if strcmp(NetStruct,'SF')
    gamma = 3;
    parameters = [gamma, k];
else
    parameters = k;
end
Anw = BuildNetwork(N, NetStruct,parameters,'gcc'); % Adjacency matrix, not weighted
n = size(Anw,1);
kreal = mean(sum(Anw));
k2 = mean(sum(Anw).^2);
kappa = k2/kreal-1;

%% initial/fixed condition
x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,holding_value);

%% simulate
toplot = 'time';
A = w*Anw;
Factor_Xeff = sum(A,1)/sum(sum(A));

pnls = 'bcd';
for it = 1:length(TimeProBioVec)
    it
    conditions.TimeProBio = TimeProBioVec(it);
    Xss = SolveOdesProbiotic(x0, M,all_or_eff, toplot,conditions,release);
            
    % save the figure in the folder 'output'
    folder = '..\..\output\FiguresSI\';
    filename = ['9',pnls(it)];
    save_pdf_min_size([folder,filename])
    close(gcf)

end

