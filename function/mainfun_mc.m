function result = mainfun_mc(R,N,f,k,d,mdl,starv)
m = length(N);
Nn = reshape(repmat(N,R,1),[],1);
Rn = R*m;
fstarv = f(starv);
%
tic
parfor i = 1:Rn
    nn = Nn(i);
    xxiid(i,:) = MCMC(nn,d,mdl,k,f);
    % disp(['finish',num2str(i)])
end
toc
for j = 1:m 
    temp = xxiid(1+(j-1)*R:j*R,:);
    rmean(j,:) = mean(temp,1);
    if flag == 1
        rmse(j,:) = mean((temp-fstarv).^2);
    else
        rvarr(j,:) = var(temp,0);
    end
end
if flag == 1
    for k = 1:d
        result(1:m,1,k) = rmean(1:m,k);
        result(m+1:2*m,1,k) = rmse(1:m,k);
    end
else
    for k = 1:d
        result(1:m,1,k) = rmean(1:m,k);
        result(m+1:2*m,1,k) = rvarr(1:m,k);
    end
end
end