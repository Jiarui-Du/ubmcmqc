%%
clc
clear

% This file implements the loading of corresponding dataset and model

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

switch modelname
    case 'Linear_boston'
        % (506*14) boston house data from R package MASS
        data = load('./data/bostondata.txt');
        Y = log(data(:,end));
        X = data(:,1:end-1);
        X(:,5) = (data(:,5)).^2;
        X(:,6) = (data(:,6)).^2;
        X(:,8) = log(data(:,8));
        X(:,9) = log(data(:,9));
        X(:,13) = log(data(:,13));
        m = length(X);
        X = [ones(m,1),X];
        p = size(X,2);
        d = p+1;
        b0 = zeros(p,1);
        B0 = 100*eye(p);
        n0 = 5;
        s0 = 0.01;
        mdl = model_linear(d,Y,X,b0,B0,n0,s0);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","X","d","p","mdl");

    case 'Linear_california'
        % (20640*9) california data from python sklearn.datasets fetch_california_housing
        Data = load('./data/california_housing_data.txt');
        data = Data;
        households = round(Data(:,5)./Data(:,6)); % households
        data(:,3) = round(Data(:,3).*households); % total_rooms
        data(:,4) = round(Data(:,4).*households); % total_bedrooms
        colNames = ["median_income","housing_median_age","total_rooms","total_bedrooms",...
            "population","avePopulation","latitude","longitude"];
        Y = log(data(:,end).*100000);
        y = log(data(:,end));
        m = length(y);
        X(:,1) = ones(m,1);
        X(:,2:4) = [data(:,1),data(:,1).^2,data(:,1).^3];% Income
        X(:,5) = log(data(:,2));
        X(:,6) = log(data(:,3)./data(:,5));
        X(:,7) = log(data(:,4)./data(:,5));
        X(:,8) = log(data(:,6));
        X(:,9) = log(households);
        p = size(X,2);
        d = p+1;
        b0 = zeros(p,1);
        B0 = 100*eye(p);
        n0 = 5;
        s0 = 0.01;
        mdl = model_linear(d,Y,X,b0,B0,n0,s0);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","X","d","p","mdl");

    case 'Probit_Vaso'
        % 39*3(Vaso) from finney
        v = [3.7,3.5,1.25,0.75,0.8,0.7,0.6,1.1,0.9,0.9,0.8,0.55,0.6,1.4,0.75,2.3,3.2,0.85,...
            1.7,1.8,0.4,0.95,1.35,1.5,1.6,0.6,1.8,0.95,1.9,1.6,2.7,2.35,1.1,1.1,1.2,0.8,0.95,0.75,1.3];
        r = [0.825,1.09,2.5,1.5,3.2,3.5,0.75,1.7,0.75,0.45,0.57,2.75,3.0,2.33,3.75,1.64,...
            1.6,1.415,1.06,1.8,2.0,1.36,1.35,1.36,1.78,1.5,1.5,1.9,0.95,0.4,0.75,0.03,1.83,2.2,2.0,3.33,1.9,1.9,1.625];
        y = [1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1];
        Y = y(:);
        M = length(v);
        X = [ones(M,1),v(:),r(:)];
        p = size(X,2);
        d = M+p;
        mdl = model_probit(d,Y,X);
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","X","d","p","mdl");

    case 'Probit_Mroz'
        % 753*8(mroz) from R-package "wooldridge" data(mroz)
        data = readtable('./data/mrozdata.csv');
        data = table2array(data);
        data(isnan(data)) = 0;
        Y = data(:,1);
        n = length(Y);
        % data(:,20) nwifeinc 6 educ 19 exper 22 exper2 5 age 3 kidslt6 4 kidsge6
        X = [ones(n,1),data(:,20),data(:,6),data(:,19),data(:,22),data(:,5),data(:,3),data(:,4)];
        zX = [X(:,1),normalize(X(:,2:end))];
        p = size(X,2);
        d = p+n;
        mdl = model_probit(d,Y,X);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"Y","X","d","p","mdl");
        
    case 'logistic_pima'
        % (392*9)Machine Learning Repository. Pima indians diabetes data set, 2012d.
        load('./data/pima.mat') 
        xx = [];
        for i = 2:8
            xx(:,i-1) = X(:,i)==0; 
        end
        Ix = ~any(xx,2);
        DataX = X(Ix,:); % Eliminating outliers
        [m,~] = size(DataX);
        dataX = [ones(m,1),DataX];
        dataY = y(Ix,:);
        [m,p] = size(dataX);
        zX = [dataX(:,1),normalize(dataX(:,2:end))];
        d = m+p;
        b = zeros(p,1);
        B = 10*eye(p);
        mdl1 = model_polygamma(d,p,dataY,zX,b,B,"mix1");
        mdl2 = model_polygamma(d,p,dataY,zX,b,B,"mix2");
        mdl = model_polygamma(d,p,dataY,zX,b,B,"direct");
        mdlIID = model_polygammaIID(d,p,dataY,zX,b,B);
        %
        filename = ['./data/',modelname,'.mat'];
        save(filename,"dataY","zX","d","p","mdl","mdlIID","mdl1","mdl2");

    case 'logistic_german'
        % (1000*49) Machine Learning Repository. Statlog (german credit data) data set, 2012b.
        DataX = load('./data/germandataX.txt');
        DataY = load('./data/germandataY.txt');
        dataX = DataX;
        dataY = DataY;
        [m,p] = size(dataX);
        zX = [dataX(:,1),normalize(dataX(:,2:end))];
        d = m+p;
        b = zeros(p,1);
        B = 10*eye(p);
        mdl = model_polygamma(d,p,dataY,zX,b,B,"direct");
        mdlIID = model_polygammaIID(d,p,dataY,zX,b,B);
        filename = ['./data/',modelname,'.mat'];
        save(filename,"dataY","zX","d","p","mdl","mdlIID");
end

