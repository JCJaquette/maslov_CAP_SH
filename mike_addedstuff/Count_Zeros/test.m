function out = test(func,xbar,gmin,gmax)

    maxr = max(xbar - gmin, gmax - xbar);

    intOfInterest = infsup(0,maxr);
    
    out = func(intOfInterest);

end

