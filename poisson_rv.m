function x=poisson_rv(lambda, N)

% % want mult*lambda >>50 -> mult=50/lambda
% % want 1/mult < 1 -> mult > 1
% x=zeros(size(lambda));
% mult=max(max(10,50./lambda));
% if any(lambda>0)
%     x(lambda>0)=round(sum(rand(ceil(mult*lambda(lambda>0)),N)<(1/mult)));
% end
% 
%  

x=gammaincinv(randn(N,1), lambda);