%%
clc
clear

% This file implements the run time of different approach of update function 
% for logistic model with the Pima dataset.

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model,
% then run 'pilot_run_tau' to simulation \tau 
% and determine the value of k.

% add function
addpath('./function')

%
modelname = 'logistic_pima';

% load data and model
filename = ['./data/',modelname,'.mat'];
load(filename)

% load k and tau
filename = ['./result/',modelname,'_tau.mat'];
load(filename)

%
f = @(x) x;
index = [10,13,16];
R = 100;
elapsedTime = zeros(R,4,length(index));
for j = 1:length(index)
    for i = 1:R
        n = index(j);
        N = 2^n;
        
        tstart0 = tic;
        ubMCMC_mean(N,d,mdlIID,k,f);
        elapsedTime(i,1,j) = toc(tstart0);

        tstart1 = tic;
        ubMCQMC_mean(N,d,"HaraseF2",mdl1,n,k,f);
        elapsedTime(i,2,j) = toc(tstart1);

        tstart2 = tic;
        ubMCQMC_mean(N,2*d-p,"HaraseF2",mdl2,n,k,f);
        elapsedTime(i,3,j) = toc(tstart2);

        tstart3 = tic;
        ubMCQMC_mean(N,d,"HaraseF2",mdl,n,k,f);
        elapsedTime(i,4,j) = toc(tstart3);

    end
    disp([num2str(n),'-finish!'])
end
meanTime = permute(mean(elapsedTime),[3,2,1]);
%
filename = ['./result/',modelname,'-diffmethod-time.mat'];
save(filename,"elapsedTime","meanTime")

