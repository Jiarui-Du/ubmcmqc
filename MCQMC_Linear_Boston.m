% Linear_boston
clc
clear

% load data and model
modelname = 'Linear_boston';
filename = ['./data/',modelname,'.mat'];
load(filename)

% load ubMCQMC result
filename = ['./result/',modelname,'.mat'];
load(filename)

%
sequence = ["iid","HaraseF2"];
n = 10:16;
R = 100;
N = 2.^n; 
f = @(x) x;
k = 5000;
starv = reshape(results(7,2,:),1,d); 
resultbiased = mainfun_qmc(n,R,N,sequence,f,k,d,mdl,starv);

filename = ['./result/biased_',modelname,'.mat'];
save(filename,'resultbiased','N','R','sequence')

%
ld = size(results,1);
m = ld/3;
result_ubmse = results(1*m+1:2*m,:,:)./R;
result_bmse = resultbiased(1*m+1:2*m,:,:)./R;
result_mse = [result_bmse,result_ubmse];
smv = log2(sqrt(result_mse));
result_tvar = sum(result_mse(:,:,1:p),3);
result_trmse = sqrt(result_tvar);
tsmv = log2(result_trmse);

% figure
style = ["o","p","x","s","d","*"];
figure
plot(n,tsmv(:,1),'b-','marker',style(1),'linestyle','-','LineWidth',2)
hold on
plot(n,tsmv(:,2),'g-','marker',style(3),'linestyle','-','LineWidth',2)
hold on
plot(n,tsmv(:,4),'r-','marker',style(4),'linestyle','-','LineWidth',2)
hold on
tn1 = n(1) + 2*(tsmv(1,1));
tn3 = n(1) + tsmv(1,4);
plot(n,-1/2*(n-tn1),'color','k','marker','none','linestyle','--','LineWidth',1.5)
hold on
plot(n,-1*(n-tn3),'color','k','marker','none','linestyle','--','LineWidth',1.5)
hold on
xlabel('log_2(N)');  
ylabel('log_2(RMSE)');
ylim([tsmv(end,4)-3 tsmv(1,1)+1])
legend(["MCMC","MCQMC-H","ubMCQMC-H","Slope:-1/2","Slope:-1"],...
    'Fontsize',16,'Location','southwest','NumColumns',2)
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))