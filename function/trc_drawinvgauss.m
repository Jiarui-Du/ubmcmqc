function xx = trc_drawinvgauss(u,t,mu,lambda,err)
% reference: Giner, G., & Smyth, G. K. (2016). statmod: probability calculations
% for the inverse Gaussian distribution. arXiv preprint arXiv:1603.06687.
n = length(mu);
Fa = 0;
if t == inf
    Fb = 1;
else
    Fb = cdfinvgauss(t,mu,lambda);
end
uab = Fa+(Fb-Fa).*u;
ml = 3*mu/(2*lambda);
x0 = mu.*((1+ml.^2).^(1/2)-ml); % mode
Im = mu >= 1e3;
x0(Im) = lambda/3;
% Ilt = uab <= 1e-4;
% x0(Ilt) = mu(Ilt)*lambda./(norminv(u(Ilt)).^2);
% Irt = uab >= 1-(1e-4);
% x0(Irt) = gaminv(u(Irt),lambda./mu(Irt),mu(Irt).^2/lambda);
xx = x0;
Isp = false(n,1);
e = zeros(n,1);
f = ones(n,1);
F = @(x) cdfinvgauss(x,mu(~Isp),lambda);
while(~all(Isp))
e(~Isp) = abs(F(xx(~Isp))-uab(~Isp));
Isp = e <= err;
F = @(x) cdfinvgauss(x,mu(~Isp),lambda);
f(~Isp) = exp(1/2*(log(lambda)-log(2*pi*xx(~Isp).^3))-lambda*((xx(~Isp)-mu(~Isp)).^2)./(2*mu(~Isp).^2.*xx(~Isp)));
xx(~Isp) = xx(~Isp)-(F(xx(~Isp))-uab(~Isp))./f(~Isp);
end
xx = max(xx,eps);
end

%%
function y = cdfinvgauss(x,mu,lambda)
z = 1./mu;
b = sqrt(lambda./x).*(x.*z - 1);
a = -sqrt(lambda./x).*(x.*z + 1);
y = normcdf(b) + exp(2*lambda*z+log(normcdf(a)));
end