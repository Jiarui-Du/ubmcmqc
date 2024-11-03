function resultcell = lossefficient_mc(K1,R,k1,vinf,d,p,mdl,varargin)
f = @(x) x;

if nargin > 7
    mdlIID = varargin{1};
    m1 = length(K1);
    Nn1 = reshape(repmat(K1,R,1),[],1);
    Rn1 = R*m1;
    xxiid = zeros(Rn1,d);
    tauiid = zeros(Rn1,1);

    %
    tic
    parfor i = 1:Rn1
        nn = Nn1(i);
        [xxiid(i,:),tauiid(i)] = ubMCMC_mean(nn,d,mdlIID,k1,f);
        % disp(['finish',num2str(i)])
    end
    toc
else
    m1 = length(K1);
    n1 = ones(m1,1);
    sequence = "iid";
    Ns = length(sequence);
    Rs = repmat(reshape(repmat(sequence,R,1),[],1),m1,1);
    Rr = R*Ns;
    Nn = reshape(repmat(K1,Rr,1),[],1);
    ns = reshape(repmat(n1,Rr,1),[],1);
    Rn1 = Rr*m1;

    xxiid = zeros(Rn1,d);
    tauiid = zeros(Rn1,1);
    %
    tic
    parfor i = 1:Rn1
        nn = Nn(i);
        nns = ns(i);
        seq = Rs(i);
        [xxiid(i,:),tauiid(i)] = ubMCQMC_mean(nn,d,seq,mdl,nns,k1,f);
        % disp(['finish',num2str(i)])
    end
    toc
end
%
rvarr = [];
for j = 1:m1
    temp = xxiid(1+(j-1)*R:j*R,:);
    rvarr(j,:) = var(temp,0)./R;
end
rtvar = sum(rvarr(:,1:p),2);
%
all_tauiid = [];
for j = 1:m1
    all_tauiid(j,:) = tauiid(1+(j-1)*R:j*R);
end
costiid = [];
ineffiid = [];
for j = 1:m1
    at = all_tauiid(j,:);
    mt = mean(at);
    costiid(j) = 2*(mt-1)+mean(max(1,K1(j)+k1-at));
    format short g
    ineffiid(j) = costiid(j)*rtvar(j)/vinf;
end
%

resultcell = cell(m1,5);
resultcell(1:m1,1) = num2cell(k1);
resultcell(1:m1,2) = num2cell(K1(:));
resultcell(1:m1,3) = num2cell(costiid(:));
resultcell(1:m1,4) = num2cell(rtvar(:));
resultcell(1:m1,5) = num2cell(ineffiid(:));