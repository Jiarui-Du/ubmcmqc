function result = cal_sample(xx,Ns,m,R,starv)
if isempty(starv)
    flag = 0;
else
    flag = 1;
end
% method, variable, sample size
for j = 1:m 
    for i = 1:Ns 
        ij = (j-1)*Ns*R;
        temp = xx(1+(i-1)*R+ij:i*R+ij,:);
        rmean(i,:,j) = mean(temp,1);
        if flag == 1
            rmse(i,:,j) = mean((temp-starv).^2,1);
        else
            rvar(i,:,j) = var(temp,0);
        end
    end
end
% sample size, method, variable
dd = size(rmean,2);
for i = 1:m
    for j = 1:Ns
        for k = 1:dd
            result_mean(i,j,k) = rmean(j,k,i);
            if flag == 1
                result_mse(i,j,k) = rmse(j,k,i);
            else
                result_var(i,j,k) = rvar(j,k,i);
            end
        end
    end
end
for k = 1:dd
    for j = 1:Ns
        if flag == 1
            result_msefactor(:,j,k) = result_mse(:,1,k)./result_mse(:,j,k);
        else
            result_factor(:,j,k) = result_var(:,1,k)./result_var(:,j,k);
        end
    end
end
if flag == 1
    result = [result_mean;result_mse;result_msefactor];
    % result.tmse = sum(result_mse,3);
else
    result = [result_mean;result_var;result_factor];
    % result.tvar = sum(result_var,3);
end
end
