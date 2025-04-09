function [ab] = chebstar2fft(a,b)

    check = 0;
    if length(a(:,1)) ~= 1 
        a = a';
        b = b';
        check = 1;
    end

    chebLa = length(a);
    chebLb = length(b);

    a = [a(1),a(2:end),a(end:-1:2)];
    b = [b(1),b(2:end),b(end:-1:2)];
    exp = ceil(log2(length(a) + length(b)));
    n = 2^exp;

    a = [a(1:chebLa),zeros(1,n - length(a)),a(chebLa+1:end)];
    b = [b(1:chebLb),zeros(1,n - length(b)),b(chebLb+1:end)];

    ab = fft(ifft(a).*ifft(b));
    ab = ab(1:n/2)*n;%ifft(a) scales the coeffs by 1/n
    ab = real(ab);

    if check == 1
        ab = ab';
    end

    norm(ab - chebstar2(a,b,n/2)) 
    % check makes it so that this fn can do chebstar with a row or column but
    % careful chebstar2 only puts out rows
    % more tests:
    % a = rand(1,ceil(100*rand(1,1)));
    % b = rand(1,length(a));
    % chebstar2fft(a,b);

end

