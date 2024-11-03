function rmean = MCMC(N,d,model,k,f)
xt = zeros(N+k-1,d);
x = model.initsample();
xt(1,:) = x;
for t = 2:N+k-1
    x = model.proposal_sample(x);
    xt(t,:) = x;
end
rmean = mean(f(xt(k:end,:)),1);
end