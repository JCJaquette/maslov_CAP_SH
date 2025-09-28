function plots=plot_coeff(coeff,order)
% this function plots the log of the norms of the coefficients against
% their order. This is done to check for decay. 
% Inputs: coeff - corresponds to the coefficients computed for either the
%                 stable or unstable manifold
%         order - order of the parameterization 
% Outputs: plots = 0
%          also generates a 2D plot
    normorder=zeros(order+1,order+1); % Not sure what this does
    plotpoints=[];
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            for j = suborder-i
                vec=zeros(1,4);
                for k=1:4
                    vec(k)=coeff(i+1,j+1,k);
                end
                
                normpoint=norm(vec);
                normorder(i+1,j+1)=normpoint;
                point=[suborder;normpoint];
                plotpoints=[plotpoints,point];
            end
        end
        suborder=suborder+1;
    end
    
  plot(plotpoints(1,:), log(plotpoints(2,:)),'o');
  xlabel('Order of the coefficients')
  ylabel('Log norm of the coefficients')
  plots=0;
end