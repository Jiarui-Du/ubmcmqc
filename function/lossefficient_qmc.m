function resultcell = lossefficient_qmc(K2,R,k2,vinf,d,p,mdl,c)
f = @(x) x;
%
m2 = length(K2);
n2 = ones(m2,1);
sequence = "sobol-Liao";
Ns = length(sequence);
Rs = repmat(reshape(repmat(sequence,R,1),[],1),m2,1);
Rr = R*Ns;
Nn = reshape(repmat(K2,Rr,1),[],1);
ns = reshape(repmat(n2,Rr,1),[],1);
Rn = Rr*m2;
%
xx = zeros(Rn,d);
tau = zeros(Rn,1);
tic
parfor i = 1:Rn
    nn = Nn(i);
    nns = ns(i);
    seq = Rs(i);
    [xx(i,:),tau(i)] = ubMCQMC_mean(nn,d,seq,mdl,nns,k2,f);
    % disp(['finish',num2str(i)])
end
%
toc
%
result_var = [];
for j = 1:m2 %样本
    temp = xx(1+(j-1)*R:j*R,:);
    result_var(j,:) = var(temp,0)./R;
end
result_tvar = sum(result_var(:,1:p),2);
%
all_tau = [];
for j = 1:m2
    all_tau(j,:) = tau(1+(j-1)*R:j*R);
end
cost = [];
ineff = [];
for j = 1:m2
    at = all_tau(j,:);
    mt = mean(at);
    cost(j) = c*(2*(mt-1)+mean(max(1,K2(j)+k2-at)));
    format short g
    ineff(j) = cost(j)*result_tvar(j)/vinf;
end

resultcell = cell(m2,5);
resultcell(1:m2,1) = num2cell(k2);
resultcell(1:m2,2) = num2cell(K2(:));
resultcell(1:m2,3) = num2cell(cost(:));
resultcell(1:m2,4) = num2cell(result_tvar(:));
resultcell(1:m2,5) = num2cell(ineff(:));