% we format a, b as N+1 square matrices
% gives the (mn)th element of the star hat cauchy product
function prod = tripstarhat(a,b,c,m,n)
    zero = 0 *a(1,1);
    total_sum=zero;
    for j=0:m
        for k=0:n
            if j == 0 && k == 0
                continue
            end
            little_sum = zero;
            for l=0:j
                for q=0:k
                    if j==m && l ==m && q == n && k == n
                        % do nothing 
                    elseif j==m && k ==n && l==0 && q == 0
                        % do nothing 
                    else
                        little_sum=little_sum+b(j-l+1,k-q+1)*c(l+1,q+1);
                    end
                end
            end
            total_sum=total_sum+a(m-j+1,n-k+1)*little_sum;
        end
    end
    prod=total_sum;
end
