clear; clc; format compact; % close all
global A Anw NameOfModel

addpath('..\functions');

% model of dynamics
NameOfModel =  'MM'; %'Eco'; % 'Ising-Sch' ; % 'Simple'; %  'SIS';% 'MAK';%  'PD';%

M = KindOfDynamics( NameOfModel );
conditions.type = 'IC';
conditions.free_value = 50;
N = 1e3; 
k = 4; % average degree
w = 0.8; % links weight

% method of attacking the structure
method  = 'weights'; % options: 'links'; % 'nodes';  %  'weights';  %

colors = [0 170 170;255 0 0]/255;
figure;ax=axes; hold on;
h = gobjects(0); % returns an empty graphics object array.

% build the network
NetStruct = 'ER';
parameters = k;
Anw = BuildNetwork(N, NetStruct,parameters,'gcc'); % Adjacency matrix, not weighted
n = size(Anw,1);
kreal = mean(sum(Anw)); % measured degree
A = w*Anw; % weighted matrix
G = graph(A);


switch method
    
    % links removal
    case 'links'
        
        links = table2array(G.Edges);
        LinksList = sub2ind(size(A),links(:,1),links(:,2));
        LinksListT = sub2ind(size(A),links(:,2),links(:,1));
        weights = links(:,3);
        m = G.numedges;
        shuffle = randperm(m);
        LinksList = LinksList(shuffle);
        LinksListT = LinksListT(shuffle);
        idxs = round(linspace(1,m,101));
        xm = zeros(size(idxs));
        x = rand(n,1);
        
        for i=1:length(idxs)-1
            disp(i)
            A(LinksList(idxs(i):idxs(i+1)-1)) = 0;
            A(LinksListT(idxs(i):idxs(i+1)-1)) = 0;
            x = SolveOdes(x, M, 'all', 0 ,conditions,1);
            xm(i) = mean(x(GCC(A)));
        end
        plot(idxs/m,xm,'-o','LineWidth',2,'Color',colors(1,:))
        qc = idxs(find(xm<0.1,1))/m;
        
        for i=1:length(idxs)-1
            A(LinksList(idxs(i):idxs(i+1)-1)) = weights(idxs(i):idxs(i+1)-1);
            A(LinksListT(idxs(i):idxs(i+1)-1)) = weights(idxs(i):idxs(i+1)-1);
            x = SolveOdes(x, M, 'all', 0 ,conditions,1);
            xm(i) = mean(x(GCC(A)));
        end
        plot((idxs(end)-idxs)/m,xm,'-*','LineWidth',1.5,'Color',colors(2,:))
        
        % Nodes removal
    case 'nodes'
        
        idxs = round(linspace(1,n,101));
        xm = zeros(size(idxs));
        x = rand(n,1);
        
        for i=1:length(idxs)-1
            disp(i)
            A(idxs(i):idxs(i+1)-1,:) = 0;
            A(:,idxs(i):idxs(i+1)-1) = 0;
            x = SolveOdes(x, M, 'all', 0 ,conditions,1);
            xm(i) = mean(x(GCC(A)));
        end
        plot(idxs/n,xm,'-o','LineWidth',1.5,'Color',colors(1,:))
        qc = idxs(find(xm<0.1,1))/n;
        
        plot((idxs(end)-idxs)/n,0,'-*','LineWidth',1.5,'Color',colors(2,:))
        
        % reducing weights
    case 'weights'
        
        qvec = linspace(0,1,101);
        wvec = (1-qvec)*w;
        xm = zeros(size(qvec));
        x = rand(n,1);
        
        for i=1:length(qvec)
            disp(i)
            A = wvec(i)*Anw;
            x = SolveOdes(x, M, 'all', 0 ,conditions,1);
            xm(i) = mean(x);
        end
        figure; hold on
        plot(qvec,xm,'-o','LineWidth',2,'Color',colors(1,:))
        qc = qvec(find(xm<0.1,1));
        
        plot(qvec,0,'-*','LineWidth',2,'Color',colors(2,:))
        
end

%% figure properties
set(gca,'FontSize',20,'box','on','LineWidth',1.5,'YTick',0:10,'XTick',[0 0.5 1])
xlabel('\boldmath$q$','Interpreter','latex','FontWeight','bold','FontSize',30);
ylabel('\boldmath$\bar{x}$','Interpreter','latex','FontWeight','bold','FontSize',30);

ylim([-0.3,3.5])
% axis square

%% save the figure in the folder 'output'
folder = '..\..\output\Figure3\';
filename = 'b';
save_pdf_min_size([folder,filename])

