clc
clear

% This file implements the simulation of \tau

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model

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

filename = ['./data/',modelname,'.mat'];
load(filename)
R1 = 1000;

if contains(modelname,'logistic')

    parfor i = 1:R1
        [at(i),~,~] = ubMCMC(1,d,mdlIID,1);
    end

    % choose k
    atau = prctile(at,99);
    k = round(2*atau);

    % save k and tau
    filename = ['./result/',modelname,'_tau.mat'];
    save(filename,'R1','at','k')
else
    sequence = "iid";
    parfor i = 1:R1
        [at(i),~,~] = ubMCQMC(1,d,sequence,mdl,1,1);
    end

    % choose k
    atau = prctile(at,99);
    k = round(2*atau);

    % save k and tau
    filename = ['./result/',modelname,'_tau.mat'];
    save(filename,'R1','at','k')
end

% plot
figure
histogram(at)
