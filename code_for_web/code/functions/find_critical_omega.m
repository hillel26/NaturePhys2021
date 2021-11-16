function [ omega ] = find_critical_omega(x0,M,toplot)
%FIND_CRITICAL_OMEGA by bisection method
global A Anw Factor_Xeff

if toplot == 0
    down=eps;
    up=1;
    while up-down>1e-4
        omega = mean([down,up]);
        A = omega*Anw;
        Factor_Xeff = sum(A,1)/sum(sum(A));
        Xss = SolveOdes(x0,M,'eff',0);
        if Xss>0.2
            up = omega;
        else
            down=omega;
        end
    end
elseif toplot==1
    down=eps;
    up=1;
    OmegaVec = zeros(1,100);
    XssVec = zeros(size(OmegaVec));
    i=0;
    while up-down>1e-4
        i=i+1;
         omega = mean([down,up]);
        A = omega*Anw;
        Factor_Xeff = sum(A,1)/sum(sum(A));
        Xss = SolveOdes(x0,M,'eff',0);
        if Xss>0.2
            up = omega;
        else
            down=omega;
        end
        OmegaVec(i) = omega;
        XssVec(i) = Xss;
    end
    figure; hold on;    
    plot(OmegaVec(1:i),XssVec(1:i),'.b');    
    x0eff = Factor_Xeff*x0;
    plot(OmegaVec(1:i),ones(1,i)*x0eff,'.r')
    legend('Steady state','Initial condition')
    xlabel('\omega')
    ylabel('x_{eff}')
    xlim([0,3*omega])
    set(gca,'XTick',omega,'XGrid','on')    
else
    disp(' --Error--: third input should be ''1'' or ''0'' to plot or not');
end
end

