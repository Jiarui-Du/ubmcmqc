function X = randinvgauss(u, mu, lambda)
% generate random variates from the inverse Gaussian distribution by
% uniform variates

% reference:Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). 
% Generating random variates using transformations with multiple roots. 
% The American Statistician, 30(2), 88-90.

p = length(u);
n = length(mu);
V = norminv(u).^2;
W = mu.*(1+mu.*V/(2*lambda));
X = W-sqrt(W.^2-mu.^2);
IU = rand(p,1) >= mu./(mu+X);
if n == 1
    X(IU) = mu^2./X(IU);
else
    X(IU) = mu(IU).^2./X(IU);
end
end
