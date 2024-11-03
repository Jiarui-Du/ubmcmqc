function [tau,xt,yt] = ubMCMC(N,d,mdl,k)
m = N+k-1;
xt = zeros(m+1,d);
yt = zeros(m,d);
x = mdl.initsample();
y = mdl.initsample();
xt(1,:) = x;
yt(1,:) = y;
x = mdl.proposal_sample(x);
xt(2,:) = x;
t = 2;
tau = inf;
while t <= max(m,tau)
    if t < tau
        [x,y] = mdl.proposal_coupling_sample(x,y);
        if x == y
            tau = t;
        end
    else
        x = mdl.proposal_sample(x);
        y = x;
    end
    xt(t+1,:) = x;
    yt(t,:) = y;
    t = t+1;
end
end