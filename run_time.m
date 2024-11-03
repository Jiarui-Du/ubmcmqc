clc
clear

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

f = @(x) x;
R = 100;
index = [10,13,16];
elapsedTime = zeros(R,2,length(index));

for j = 1:length(index)
    for i = 1:R
        n = index(j);
        N = 2^n;

        if contains(modelname,'logistic')
            tstart1 = tic;
            ubMCMC_mean(N,d,mdlIID,k,f);
            elapsedTime(i,1,j) = toc(tstart1);
        else
            tstart1 = tic;
            ubMCQMC_mean(N,d,"iid",mdl,n,k,f);
            elapsedTime(i,1,j) = toc(tstart1);
        end

        tstart2 = tic;
        ubMCQMC_mean(N,d,"HaraseF2",mdl,n,k,f);
        elapsedTime(i,2,j) = toc(tstart2);

    end
    disp([num2str(n),'-finish!'])
end

meanTime = permute(mean(elapsedTime,1),[3,2,1]);

filename = ['./result/Time-',modelname,'.mat'];
save(filename,'elapsedTime','meanTime')

disp([modelname,' finish!'])
