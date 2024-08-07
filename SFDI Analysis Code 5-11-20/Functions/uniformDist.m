function x = uniformDist(lower,upper,N)
    % uniformDist returns N values on the open interval (lower,upper)
    
    % transform (0,1) to (lower,uppper) using rand
    x = (upper - lower)*rand(1,N) + lower;
    
end %function