function plots=plot_coeff_sum(coeff,params,marker)
% this function plots the log of the norms of the coefficients against
% their order. This is done to check for decay. 
% Inputs: coeff - corresponds to the coefficients computed for either the
%                 stable or unstable manifold
%         order - order of the parameterization 
% Outputs: plots = 0
%          also generates a 2D plot
    % normorder=zeros(order+1,order+1); % Not sure what this does

    order=params.order;

    if params.isIntval
        coeff =coeff.mid;
    end

    suborder=0;
    normorder=zeros(order+1,2);
    while suborder<order+1
        little_sum = 0;
        for i=0:suborder
            for j = suborder-i
                little_sum = little_sum + sum(abs(coeff(i+1,j+1,:)),'all');
            end
        end
        normorder(suborder+1,:)=[suborder little_sum];
        suborder=suborder+1;
    end
    
  plot(normorder(:,1), log(normorder(:,2)) / log(10),marker);
  xlabel('Order $n$ of coefficients','Interpreter','latex')
  ylabel('$\log_{10} \sum_{|\alpha| =n} \| P_\alpha \| $ ','Interpreter','latex','Rotation',0)
  plots=0;
end