function x0 = Initial_or_fixed_Conditions(n,NumExcited,ChooseExcNodes,conditions,value)
global Anw List

% switch conditions_type
%     case 'IC'
%         product = 10*n;
%         value = product / NumExcited;
%     case 'BC'
%         value = 1;
% end
x0 = zeros(n,1);

if NumExcited>n
    fprintf ( 2, 'warning!! number of wanted excited nodes larger than the system size\n' );
    NumExcited=n;
end

switch ChooseExcNodes
    case 'rand'
        shake = randperm(n);
        IdxsNotZero = shake(1:NumExcited);
    case 'low degree'
        [~, Kidxsort] = sort(sum(Anw)) ; % for high or low degree excited nodes
        IdxsNotZero = Kidxsort(1:NumExcited);
    case 'high degree'
        [~, Kidxsort] = sort(sum(Anw),'descend') ; % for high or low degree excited nodes
        IdxsNotZero = Kidxsort(1:NumExcited);
    case 'local'
        % BFS for localized excitation
        CenterNode = randi(n);
            GenerationsArray = graphtraverse(Anw,CenterNode,'directed','false','Method','Bfs');
%         GenerationsArray = BFSforGivenNode(List,CenterNode);
        if length(GenerationsArray)<NumExcited
            CenterNode = randi(n);
            GenerationsArray = graphtraverse(Anw,CenterNode,'directed','false','Method','Bfs');
        end
        IdxsNotZero = GenerationsArray(1:min(NumExcited,end));
    case '1'
        IdxsNotZero = 1;
    case 'specific' % given array of nodes
        IdxsNotZero = conditions.specific_nodes(1:NumExcited);
end

x0(IdxsNotZero) = value;

end


