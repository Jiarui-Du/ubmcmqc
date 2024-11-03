function rmean = MCQMC(N,d,sequence,model,n,k,f)
xt = zeros(N+k-1,d);
x = model.initsample();
xt(1,:) = x;
u = seqfun(sequence,N,d,n);
u = [rand(k-1,d);u];
for t = 1:N+k-1
    u1 = u(t,:);
    x = model.proposal_sample(x,u1);
    xt(t+1,:) = x;
end
rmean = mean(f(xt(k:end,:)),1);
end