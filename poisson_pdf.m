function P=poisson_pdf(x, lambda)


% originally: P=lambda^k.*exp(-lambda)/factorial(k).  This blows up!
% take the log:
% log(P) = k*log(lambda)-lambda -log(factorial(x))
% and factorial(x) = gamma(x+1)
P=exp(-lambda+x.*log(lambda)-gammaln(x+1));