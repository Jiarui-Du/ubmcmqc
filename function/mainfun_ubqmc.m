function results = mainfun_ubqmc(n,R,N,sequence,f,k,d,mdl,starv)
m = length(N);
Ns = length(sequence);
Rs = repmat(reshape(repmat(sequence,R,1),[],1),m,1);
Rr = R*Ns;
Nn = reshape(repmat(N,Rr,1),[],1);
ns = reshape(repmat(n,Rr,1),[],1);
Rn = Rr*m;
fstarv = f(starv);
%
tic
parfor i = 1:Rn
    nn = Nn(i);
    nns = ns(i);
    seq = Rs(i);
    xx(i,:) = ubMCQMC_mean(nn,d,seq,mdl,nns,k,f);
    % disp(['finish',num2str(i)])
end
results = cal_sample(xx,Ns,m,R,fstarv);
end