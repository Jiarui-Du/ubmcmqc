clc
clear

% This file implements the RMSE result and convergence plot.

% One should run 'data_model.m' at first 
% to load the corresponding dataset and model,
% then run 'pilot_run_tau' to simulation \tau 
% and determine the value of k.

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

n = 10:16;
m = length(n);
R = 100;
N = 2.^n;
f = @(x) x;
starv = [];

if contains(modelname,'logistic')
    % parallel run different N, different sequence

    resultIID = mainfun_ubmc(R,N,f,k,d,mdlIID,starv);

    sequence = ["HaraseF2","sobol-Liao"];
    results = mainfun_ubqmc(n,R,N,sequence,f,k,d,mdl,starv);

    result_mean = [resultIID(1:m,:,:),results(1:m,:,:)];
    result_var = [resultIID(1*m+1:2*m,:,:),results(1*m+1:2*m,:,:)]./R;
    for i = 1:3
        result_fac(:,i,:) =  result_var(:,1,:)./result_var(:,i,:);
    end
    % row : N from 2^10 to 2^16
    % column : sequence with iid, Harase, Liao
    
    % save result
    filename = ['./result/',modelname,'.mat'];
    save(filename,'results','resultIID','result_var','N','R')

    disp([modelname,' finish!'])
else

    % parallel run different N, different sequence
    sequence = ["iid","HaraseF2","sobol-Liao"];
    results = mainfun_ubqmc(n,R,N,sequence,f,k,d,mdl,starv);
    % row : N from 2^10 to 2^16
    % column : sequence with iid, Harase, Liao
    % row 1-m mean;m+1-2m variance; 2m+1-3m factor
    result_mean = results(1:m,:,:);
    result_var = results(1*m+1:2*m,:,:)./R;
    result_fac = results(2*m+1:3*m,:,:);

    % save result
    filename = ['./result/',modelname,'.mat'];
    save(filename,'results','result_var','N','R')

    disp([modelname,' finish!'])

end
format short g

% table and figure
% total_var
result_mse = result_var;
result_rmse = sqrt(result_var);
result_tmse = sum(result_mse(:,:,1:p),3);
mm = size(result_tmse,2);
result_tmrf = zeros(m,mm);
for i = 1:mm
    result_tmrf(:,i) = result_tmse(:,1)./result_tmse(:,i);
end
result_trmse = sqrt(result_tmse);
result_trrf = zeros(m,mm);
for i = 1:mm
    result_trrf(:,i) = result_trmse(:,1)./result_trmse(:,i);
end
index = [1,4,7];
% table for presentation
table = [result_trmse(index,1),result_trrf(index,2)]';

smv = log2(result_rmse);
tsmv = log2(result_trmse);

%% convergence rate
k1 = polyfit(n,tsmv(:,1),1);
k2 = polyfit(n,tsmv(:,2),1);
k3 = polyfit(n,tsmv(:,3),1);

% convergence plot of total rmse
style = ["o","p","x","s","d","*"];
figure
n = 10:16;
plot(n,tsmv(:,1),'b-','marker',style(1),'linestyle','-','LineWidth',2)
hold on
plot(n,tsmv(:,2),'r-','marker',style(2),'linestyle','-','LineWidth',2)
hold on
plot(n,tsmv(:,3),'k-','marker',style(3),'linestyle','-','LineWidth',2)
hold on
tn1 = n(1) + 2*(tsmv(1,1));
tn2 = n(1) + tsmv(1,2);
tn3 = n(1) + 2*tsmv(1,3);
tn4 = n(1) + tsmv(1,3);
plot(n,-1/2*(n-tn1),'color','k','marker','none','linestyle','--','LineWidth',1.5)
hold on
plot(n,-1*(n-tn2),'color','k','marker','none','linestyle','-.','LineWidth',1.5)
hold on
plot(n,-1/2*(n-tn3),'color','k','marker','none','linestyle','--','LineWidth',1.5)
hold on
plot(n,-1*(n-tn4),'color','k','marker','none','linestyle','-.','LineWidth',1.5)
xlabel('log_2(N)');
ylabel('log_2(RMSE)');
ylim([tsmv(end,2)-3.5 tsmv(1,1)+1])
legend(["ubMCMC","ubMCQMC-H","ubMCQMC-L","Slope:-1/2","Slope:-1"],...
    'Fontsize',16,'Location','southwest','NumColumns',2)
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))

% all interesting components
style = ["o","p","x","s","d","*"];
color = ["b","r","k"];
figure
ID = 1:p;
index = [1,4,7];
for i = index
    plot(ID,reshape(smv(i,1,ID),1,p),'color',color(1),'marker',style((i+2)/3),'linestyle','-','LineWidth',2)
    hold on
end
for i = index
    plot(ID,reshape(smv(i,2,ID),1,p),'color',color(2),'marker',style((i+2)/3),'linestyle','-','LineWidth',2)
    hold on
end
xlabel('Components');
ylabel('log_2(RMSE)');
xlim([1 p])
ylim([min(reshape(smv(index(3),2,ID),1,p))-12 floor(max(reshape(smv(index(1),1,ID),1,p))+2)])
legend(["ubMCMC-2^{10}","ubMCMC-2^{13}","ubMCMC-2^{16}",...
    "ubMCQMC-H-2^{10}","ubMCQMC-H-2^{13}","ubMCQMC-H-2^{16}"],'Fontsize',16,'Location','south','NumColumns',2)
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))

