clc
clear

% This file implements the impact of different k for ubMCQMC

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model,
% then run 'pilot_run_tau' to simulation \tau 

% add function
addpath('./function')

% choose model
allmodel = ["Linear_boston","Linear_california","Probit_Vaso",...
    "Probit_Mroz","logistic_pima","logistic_german"];

modelnumber = 1; % different number for different model
% 1: Linear_boston
% 2: Linear_california
% 3: Probit_Vaso
% 4: Probit_Mroz
% 5: logistic_pima
% 6: logistic_german
modelname = char(allmodel(modelnumber));

% load data and model
filename = ['./data/',modelname,'.mat'];
load(filename)

% load k and tau
filename = ['./result/',modelname,'_tau.mat'];
load(filename)
atau = prctile(at,99);
k = round(atau);

% repeat 25 times
r = 25;
sequence = ["iid","HaraseF2","sobol-Liao"];
n = 10;
R = 100;
N = 2.^n;
f = @(x) x;
starv = [];
%
K = [1,k,2*k,4*k];
result1_tvar = zeros(r,3,4);
result2_tvar = zeros(r,3,4);

for i = 1:4
    % case 1 direct using WCUD
    [~,~,result1_tvar(:,:,i)] = mainfun_ubqmc_k(n,r,r*R,N,sequence,f,K(i),d,p,mdl,starv,0);
    % case 2 using WCUD after burnin period
    [~,~,result2_tvar(:,:,i)] = mainfun_ubqmc_k(n,r,r*R,N,sequence,f,K(i),d,p,mdl,starv);
end

filename = ['./result/',modelname,'-diffk','.mat'];
save(filename,"result1_tvar","result2_tvar","r","R","k")

disp([modelname,' finish!'])

% filename = ['./result/',modelname,'-diffk','.mat'];
% load(filename)

result1_trmse = sqrt(result1_tvar);
t1mean = mean(result1_trmse);
result1_tmean = (reshape(t1mean,3,4))';
result1_tfmean = zeros(4,3);
for kk = 1:3
    result1_tfmean(:,kk) = result1_tmean(:,1)./result1_tmean(:,kk);
end

t1cv = sqrt(var(result1_trmse,0))./t1mean;
result1_tcv = (reshape(t1cv,3,4))';
result1_tfcv = zeros(4,3);
for kk = 1:3
    result1_tfcv(:,kk) = result1_tcv(:,1)./result1_tcv(:,kk);
end

table1meancv = [result1_tmean(:,1),result1_tcv(:,1)];
table1meancv = [table1meancv,result1_tfmean(:,2),result1_tcv(:,2),...
    result1_tfmean(:,3),result1_tcv(:,3)];

result2_trmse = sqrt(result2_tvar);
t2mean = mean(result2_trmse);
result2_tmean = (reshape(t2mean,3,4))';
result2_tfmean = zeros(4,3);
for kk = 1:3
    result2_tfmean(:,kk) = result2_tmean(:,1)./result2_tmean(:,kk);
end

t2cv = sqrt(var(result2_trmse,0))./t2mean;
result2_tcv = (reshape(t2cv,3,4))';
result2_tfcv = zeros(4,3);
for kk = 1:3
    result2_tfcv(:,kk) = result2_tcv(:,1)./result2_tcv(:,kk);
end

table2meancv = [result2_tmean(:,1),result2_tcv(:,1)];
table2meancv = [table2meancv,result2_tfmean(:,2),result2_tcv(:,2),...
    result2_tfmean(:,3),result2_tcv(:,3)];

table = [table1meancv;table2meancv];