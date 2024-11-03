function X = pgdrawnew(u,z)
% X = pgdrawnew() generates random variates from PG(1, z) by unifrom variates

% references: 

% Makalic, E., & Schmidt, D. F. (2016). 
% High-dimensional Bayesian regularised regression with the BayesReg package. 
% arXiv preprint arXiv:1611.06649.

% Polson, N. G., Scott, J. G., & Windle, J. (2013). 
% Bayesian inference for logistic models using Pólya–Gamma latent variables.
% Journal of the American statistical Association, 108(504), 1339-1349.

% PhD thesis of Windle (2013)

% J*(1, z)
n = length(z);
t = 2/pi;
% t = 0.64;
z = abs(z)/2;
%% computer the ratio q / (p + q) （Propobility of inverse Gaussian distribution）
% Note that the p and q here are in the opposite form of the paper Polson et al.(2013), 
% Windle convert the two symbols in the pseudocode
K = z.^2/2 + pi^2/8;
Kt = K*t;
lambda = 1;
mu = 1./z;
logp = log(pi)-log(2*K)-Kt;
logq = log(2)-z+log(cdfinvgauss(t,mu,lambda));
ratio = 1./(1+exp(logp-logq));
%% setup variables for vectorisation
X = zeros(n,1);
Isampled = false(n,1);
Sn = zeros(n,1);
% uu = zeros(n,1);
%%
while(~all(Isampled))
%% Step 1: Sample X ? g(x|z)
%% Sample exponential when F(x) > ratio
Ixe = ~Isampled & (u > ratio);
X(Ixe) = t-log(1-(u(Ixe)-ratio(Ixe))./(1-ratio(Ixe)))./K(Ixe);
%% Sampled from truncated inverse-gaussian when F(x) < ratio
Ixig = ~Isampled & (u < ratio);
err = 1e-14;
sampleX = trc_drawinvgauss(u(Ixig)./ratio(Ixig),t,mu(Ixig),lambda,err);
X(Ixig) = sampleX;
%% Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U <= Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    i = 1;
    Ix = ~Isampled;
    Sn(Ix) = an(0, X(Ix), t);
    U = rand(n,1) .* Sn;
    even = false;
    sign = -1;
    
    Idone = Isampled;
    %%
    % While some RVs have either not been accepted, or have not yet been
    % rejected ...
    while(~all(Idone))
        Ix = ~Idone;
        Sn(Ix) = Sn(Ix) + sign*an(i, X(Ix), t);
        
        %Sn = Sn + sign * a(i, X, t);
        % Accept if n is odd
        Ixac = ((U <= Sn) & ~even & ~Isampled);
        X(Ixac) = X(Ixac)/4;
        Isampled(Ixac) = true;
        Idone(Ixac) = true;
        
        % Return to step 1 if n is even
        Ixre = (U > Sn) & ~Isampled & ~Idone & even;
        X(Ixre) = X(Ixre);
        % u = ur2;
        u(Ixre) = rand(sum(Ixre),1);
        Idone(Ixre) = true;
        even = ~even;
        sign = -sign;
        i = i + 1;
    end
end
end
% 
% Function a_n(x) defined in equations (12) and (13) of
% Bayesian inference for logistic models using Polya-Gamma latent variables
% Nicholas G. Polson, James G. Scott, Jesse Windle
% arXiv:1205.0310
%
% Also found in the PhD thesis of Windle (2013) in equations
% (2.14) and (2.15), page 24

function f = an(n,x,t)
f = zeros(length(x), 1);
Ix = x <= t;
f(Ix) = log(pi) + log(n + 0.5) + (3/2)*(log(2)-log(pi)-log(x(Ix))) - 2*(n + 0.5).^2./x(Ix);
Ix = x > t;
f(Ix)  = log(pi) + log(n + 0.5) - x(Ix) * pi^2 / 2 * (n + 0.5)^2;
f = exp(f);
end

function y = cdfinvgauss(x,mu,lambda)
z = 1./mu;
b = sqrt(lambda./x).*(x.*z - 1);
a = -sqrt(lambda./x).*(x.*z + 1);
y = normcdf(b) + exp(2*lambda*z+log(normcdf(a)));
end
% 
