function result=gaussian(x,c,sigma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% result=gaussian(x,c,sigma)
% returns the value of a gaussian 
% with standard deviation sigma,
% centered at c, evaulated at x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result=zeros(size(x));
result=1/(sigma*sqrt(2*pi)).*exp(-(x-c).^2/2/sigma^2);


