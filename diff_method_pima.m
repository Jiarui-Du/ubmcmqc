%%
clc
clear

% This file implements different approach of update function 
% for logistic model with the Pima dataset.

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model,
% then run 'pilot_run_tau' to simulation \tau 
% and determine the value of k.

%
modelname = 'logistic_pima';

% load data and model
filename = ['./data/',modelname,'.mat'];
load(filename)

% load k and tau
filename = ['./result/',modelname,'_tau.mat'];
load(filename)

% IID
n = 10:16;
m = length(n);
R = 100;
N = 2.^n; 
f = @(x) x;
starv = []; 
resultIID = mainfun_ubmc(R,N,f,k,d,mdlIID,starv);

% CUD
sequence = ["HaraseF2","sobol-Liao"];
results1 = mainfun_ubqmc(n,R,N,sequence,f,k,d,mdl1,starv);
results2 = mainfun_ubqmc(n,R,N,sequence,f,k,2*d-p,mdl2,starv);
results3 = mainfun_ubqmc(n,R,N,sequence,f,k,d,mdl,starv);

% save
filename = ['./result/',modelname,'-diffmethod.mat'];
save(filename,'resultIID','results1','results2','results3','N','sequence')

disp([modelname,' finish!'])

%
resultiid_mean = resultIID(1:m,:,:);
resultiid_var = resultIID(1*m+1:2*m,:,:)./R;
result1_mean = results1(1:m,:,:);
result1_var = results1(1*m+1:2*m,:,:)./R;
result2_mean = results2(1:m,:,:);
result2_var = results2(1*m+1:2*m,:,:)./R;
result3_mean = results3(1:m,:,:);
result3_var = results3(1*m+1:2*m,:,:)./R;

result_mean = [resultiid_mean,result1_mean,result2_mean,result3_mean];
result_var = [resultiid_var,result1_var,result2_var,result3_var];
mm = 1+2*3;
result_fac = zeros(m,mm,d);
for i = 1:mm
    result_fac(:,i,:) = result_var(:,1,:)./result_var(:,i,:);
end
%
format short g
result_tvar = sum(result_var(:,:,1:p),3);
resultall_tfac = [];
result_tfac = zeros(m,mm);
for i = 1:mm
    result_tfac(:,i) = result_tvar(:,1)./result_tvar(:,i);
end
result_trmse = sqrt(result_tvar);
result_tfrmse = zeros(m,mm);
for i = 1:mm
    result_tfrmse(:,i) = result_trmse(:,1)./result_trmse(:,i);
end
index1 = [2,4,6];
table = [result_trmse(:,1),result_tfrmse(:,index1)];
index = [1,4,7];
table2 = table(index,:);