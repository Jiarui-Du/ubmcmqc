function X = pgdraw_cudiidz(u,z)
% X = pgdraw_cudiidz() generates random variates from PG(1, z) by unifrom variates

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
%% computer the ratio p / (p + q)
K = z.^2/2 + pi^2/8;
Kt = K*t;
lambda = 1;
mu = 1./z;
logp = log(pi)-log(2*K)-Kt;
logq = log(2)-z+log(cdfinvgauss(t,mu,lambda));
ratio = 1./(1+exp(logq-logp));
%% setup variables for vectorisation
X = zeros(n,1);
Isampled = false(n,1);
Sn = zeros(n,1);
uu = zeros(n,1);
%%
while(~all(Isampled))
%% Step 1: Sample X ? g(x|z)
uu(~Isampled) = rand(sum(~Isampled),1);
%% Sample exponential, as required
Ixe = ~Isampled & (uu < ratio);
X(Ixe) = t-log(1-u(Ixe))./K(Ixe);
%% Sampled from truncated inverse-gaussian, as required
Ixig = ~Isampled & ~(uu < ratio);
% err = 1e-14;
% sampleX = drawg(u,t,mu,lambda,n,Ixig,err);
% sampleX = trc_drawinvgauss(u(Ixig),t,mu(Ixig),lambda,err);
sampleX = truncinvgrng_vec(z(Ixig),t,u(Ixig));
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

% %%
function y = cdfinvgauss(x,mu,lambda)
z = 1./mu;
b = sqrt(lambda./x).*(x.*z - 1);
a = -sqrt(lambda./x).*(x.*z + 1);
y = normcdf(b) + exp(2*lambda*z+log(normcdf(a)));
end
% 
function X = truncinvgrng_vec(z,t,u)
mu = 1./z;
n = length(mu);
X = zeros(n, 1);
%% Rejection sampler based on truncated gamma
Ix = mu > t;
nx = sum(Ix);
Z = zeros(nx, 1);
Isampled = false(nx, 1);
zz = z(Ix);
while (~all(Isampled))
    uu = rand(nx, 1);
    % Z(~Isampled) = 1.0 ./ truncgamma_vec(sum(~Isampled), 1/t);
    Z(~Isampled) = 1.0 ./ truncgamma_inv(u(~Isampled),1/2,1/2,1/t,inf);
    
    Isampled(uu < exp(-zz.^2/2 .* Z)) = true;
end
X(Ix) = Z;
%% Direct rejection sampler
Ix = ~Ix;
Isampled = false(sum(Ix), 1);
Z = zeros(sum(Ix), 1);
m = mu(Ix);
while (~all(Isampled))
    % Z(~Isampled) = randinvg(m(~Isampled), 1);
    Z(~Isampled) = randinvgauss(u(~Isampled),m(~Isampled),1);
    Isampled(Z < t) = true;
end
X(Ix) = Z;
end

function X  = truncgamma_inv(u,a,b,ia,ib)
Fa = gamcdf(ia,a,1/b);
Fb = gamcdf(ib,a,1/b);
uab = Fa+(Fb-Fa)*u;
X = gaminv(uab,a,1/b);
end

% Sample truncated gamma random variates
%
%   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities 
%   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
% function X  = truncgamma_vec(n, c)
% X = zeros(n, 1);
% Isampled = false(n, 1);
% gX = zeros(n, 1);
% while (~all(Isampled))
%     X(~Isampled) = c + exprnd_fast(ones(sum(~Isampled),1))*2;
%     gX(~Isampled) = sqrt(2/pi) ./ sqrt(X(~Isampled));%    exp(-(b-g)*X(~Isampled))/K./sqrt(X(~Isampled));
% 
%     Isampled(~Isampled) = (rand(sum(~Isampled),1) <= gX(~Isampled));
% end
% end
% 
% function r = exprnd_fast(mu)
% % length of mu
% n = length(mu);
% 
% % Generate uniform random values, and apply the exponential inverse CDF.
% r = -mu .* log(rand(n, 1));
% 
% end