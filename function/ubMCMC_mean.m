function [rmean,tau] = ubMCMC_mean(N,d,mdl,k,f)

[tau,xt,yt] = ubMCMC(N,d,mdl,k);

dd = length(f(xt(1,:)));
hkm = zeros(N,dd);
for l = 1:N
    tt = k-1+l;
    if tt >= tau
        hkm(l,:) = f(xt(tt+1,:));
    else
        hkm(l,:) = sum(f(xt(tt+1:tau,:)))-sum(f(yt(tt+1:tau-1,:)));
    end
end
rmean = mean(hkm,1);
end