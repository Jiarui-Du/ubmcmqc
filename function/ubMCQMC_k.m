function rmean = ubMCQMC_k(N,d,sequence,mdl,n,k,f)
m = N;
md = mdl.dimension;
xt = zeros(m+1,md);
yt = zeros(m,md);
x = mdl.initsample();
y = mdl.initsample();
xt(1,:) = x;
yt(1,:) = y;

u = seqfun(sequence,N,d,n);
% u = [rand(kk-1,d);u];

% tic
u1 = u(1,:);
x = mdl.proposal_sample(x,u1);
xt(2,:) = x;
t = 2;
tau = inf;
while t <= max(m,tau)
    if t <= m
        u1 = u(t,:);
    else
        u1 = rand(1,d);
    end
    if t < tau
        [x,y] = mdl.proposal_coupling_sample(x,y,u1);
        if x == y
            tau = t;
        end
    else
        x = mdl.proposal_sample(x,u1);
        y = x;
    end
    xt(t+1,:) = x;
    yt(t,:) = y;
    t = t+1;
end
dd = length(f(xt(1,:)));
M = m-k+1;
hkm = zeros(M,dd);
for l = 1:M
    tt = k-1+l;
    if tt >= tau
        hkm(l,:) = f(xt(tt+1,:));
    else
        hkm(l,:) = sum(f(xt(tt+1:tau,:)))-sum(f(yt(tt+1:tau-1,:)));
    end
end
rmean = mean(hkm,1);
end