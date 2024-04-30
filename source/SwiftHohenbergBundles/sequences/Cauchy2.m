% a, b: structs with fields 
    % int order: max order 
    % array coeff: vector of coeffs
    % int index_lb: lowest order 
    % int dim: dimension of coeffs (i.e. dim of multi-index) 
% This paper looks promising for FFT implementation -- http://www1.cs.columbia.edu/~stratos/research/fft.pdf
function coeff = Cauchy2(a,b,order)
    % assert dims are the same 
    % assert order <= min(a.order + b.index_lb, b.order + a.index_lb) 

    index_lb = a.index_lb + b.index_lb; 
   % coeff = zeros(order + 1 - index_lb,  order + 1 - index_lb); 

   coeff = zeros(order + 1, order + 1)

    % Assuming a starts with index -2, pad b out with zeros so that it also starts
    % with index -2. Then we can use the same indexing for a and b

    if a.index_lb > b.index_lb 
        temp = a; 
        a = b; 
        b = temp; 
    end 
    
    % assumes the difference between index lb of a and b is 2 
    zero_col = zeros(size(a.coeff,1), 2); 
    zero_row = zeros(2, size(a.coeff,2) + 2);
    b_half_pad = [zero_col, b.coeff]; 
    b_pad = [zero_row; b_half_pad]; 

    a_shift = -a.index_lb + 1;
    b_shift = -b.index_lb + 1;
    for m = index_lb:1:order
        for n = index_lb:1:order - m 
            % disp("mn")
            % disp(m)
            % disp(n)
            pmn = 0;
            for i = index_lb:1:m 
                for j = index_lb:1:n 
    
    
                        pmn = pmn ... 
                        + a.coeff(i+a_shift, j+a_shift)*b.coeff(m-i+b_shift, n-j+b_shift);
                end
            end
              coeff(m - index_lb + 1, n - index_lb + 1) = pmn;

            end
        end 


    end


