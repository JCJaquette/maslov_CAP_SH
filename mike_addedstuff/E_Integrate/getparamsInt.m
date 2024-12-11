function params = getparamsInt(n)

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;

    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;

    return

end

    error('Make sure n=1, 2, or 3')

end



