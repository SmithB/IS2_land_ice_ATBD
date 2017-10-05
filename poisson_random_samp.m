function N=poisson_random_samp(lambda)

if lambda==0
    N=0; return
end

if lambda > 5
    N=round(poisson_rv(lambda, 1));
else
    p=1; k=0; 
    L=exp(-lambda);
    while p>L
        p=p*rand(1);
        k=k+1;
    end
    N=k-1;
end
    
    
