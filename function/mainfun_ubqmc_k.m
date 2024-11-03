function [rmean,rvar,tvar] = mainfun_ubqmc_k(n,r,R,N,sequence,f,k,d,p,mdl,starv,varargin)

m = length(N);
Ns = length(sequence);
Rs = repmat(reshape(repmat(sequence,R,1),[],1),m,1);
Rr = R*Ns;
Nn = reshape(repmat(N,Rr,1),[],1);
ns = reshape(repmat(n,Rr,1),[],1);
Rn = Rr*m;
fstarv = f(starv);
%

if nargin > 11
    % case 1 direct using WCUD
    parfor i = 1:Rn
        nn = Nn(i);
        nns = ns(i);
        seq = Rs(i);
        xx(i,:) = ubMCQMC_k(nn,d,seq,mdl,nns,k,f);
    end
else
    % case 2 using WCUD after burnin period
    parfor i = 1:Rn
    nn = Nn(i);
    nns = ns(i);
    seq = Rs(i);
    xx(i,:) = ubMCQMC_mean(nn,d,seq,mdl,nns,k,f);
    end
end

R1 = R/r;
for i = 1:Ns 
    for j = 1:r
        ij = (i-1)*R;
        temp = xx(1+(j-1)*R1+ij:j*R1+ij,:);
        rmean(j,:,i) = mean(temp,1);
        rvar(j,:,i) = var(temp,0)./R1;
        tvar(j,i) = sum(rvar(j,1:p,i));
    end
end
end