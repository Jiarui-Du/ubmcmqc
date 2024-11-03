clc
clear

% This file implements the loss of efficiency of different method

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model,
% then run 'pilot_run_tau' to simulation \tau 
% and determine the value of k.

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

% run long chain MCMC
R = 1000;
kk = 10^3;
M = 5*10^5;
N = M-kk+1;
f = @(x) x;
meaninf = zeros(R,d);

% filename = ['./result/',modelname,'_infMCMC.mat'];
% load(filename)

if contains(modelname,'logistic')
    parfor i = 1:R
        meaninf(i,:) = MCMC(N,d,mdlIID,kk,f);
    end
else
    parfor i = 1:R
        meaninf(i,:) = MCQMC(N,d,"iid",mdl,1,kk,f);
    end
end

filename = ['./result/',modelname,'_infMCMC.mat'];
save(filename,'meaninf','N','M','kk')

% compute
rvar = sum(var(meaninf(:,1:p),0)./R);
vinf = rvar*N;

switch modelname
    case 'Linear_boston'
        K1 = [80,120,160,200];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl);
        %
        K2 = [5,6,7,8];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,1);

    case 'Linear_california'
        K1 = [60,90,120,150,180,210];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl);
        %
        K2 = [4,6,7,8,9];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,1);

    case 'Probit_Vaso'
        %
        K1 = [738,820,902,984];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl);
        %
        K2 = [2^4,20,24,28,2^5];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,1);

    case 'Probit_Mroz'
        %
        K1 = [500,550,600,650,700];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl);
        %
        K2 = [15,16,17,18,20];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,1);

    case 'logistic_pima'
        %
        K1 = [320,384,448,512,576,640];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl,mdlIID);
        %
        K2 = [2^4,20,22,24,26,28,2^5];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,2);

    case 'logistic_german'
        %
        K1 = [1300,1950,2600,3250,3900];
        resultcell_iid = lossefficient_mc(K1,R,k,vinf,d,p,mdl,mdlIID);
        %
        K2 = [60,2^6,68,72,76,80];
        resultcell_cud = lossefficient_qmc(K2,R,k,vinf,d,p,mdl,2);
end

resultcell1 = {"k","N","cost","Var","Ineff/VMC"};
resultcell = [resultcell1;resultcell_iid;resultcell_cud];
%
filename = ['./result/Ineff-',modelname,'.mat'];
save(filename,'resultcell','R')

disp([modelname,' finish!'])