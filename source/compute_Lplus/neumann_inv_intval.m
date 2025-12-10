function [outputArg1,outputArg2] = neumann_inv_intval(A)

    A_n = mid(A);
    A_r = A - mid(A);

    A_x = -A_n\A_r;



end