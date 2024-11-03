classdef model_polygammaIID
    % 使用单个变量
    % MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        para_dimension
        dataY
        dataX
        datab
        dataB
        init_sample
    end

    methods
        function obj = model_polygammaIID(dimension,para_dimension,Y,X,b,B,varargin)
            obj.ModelName  = 'polygammaGibbs';
            obj.dimension = dimension;
            obj.para_dimension = para_dimension;
            obj.dataY = Y;
            obj.dataX = X;
            obj.dataB = B;
            obj.datab = b;
            if nargin > 6
                obj.init_sample = sample0;
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function x0 = initsample(obj)
            d = obj.dimension;
            X = obj.dataX;
            B = obj.dataB;
            b = obj.datab;
            beta = mvnrnd(b,B);
            xb = X*beta(:);
            W = pgdraw(abs(xb));
            x0 = reshape([beta(:);W(:)],1,d);
        end

        %% proposal sampler
        function y = proposal_sample(obj,x)
            d = obj.dimension;
            p = obj.para_dimension;
            Y = obj.dataY;
            X = obj.dataX;
            B = obj.dataB;
            b = obj.datab;

            W = x(p+1:end);
            Sigma = inv(X'*diag(W)*X+inv(B));
            mu = Sigma*(X'*(Y-1/2)+inv(B)*b);
            beta = mvnrnd(mu,Sigma);

            Z = abs(X*beta(:));
            W = pgdraw(Z);
            y = reshape([beta(:);W(:)],1,d);
        end

        %% proposal pdf
        function [xk,yk] = proposal_coupling_sample(obj,x,y)
            d = obj.dimension;
            xk = zeros(1,d);
            yk = zeros(1,d);
            pd = obj.para_dimension;
            Y = obj.dataY;
            X = obj.dataX;
            B = obj.dataB;
            b = obj.datab;
            n = length(Y);

            Wx = x(pd+1:end);
            Sigmap = inv(X'*diag(Wx)*X+inv(B));
            mup = Sigmap*(X'*(Y-1/2)+inv(B)*b);
            betax = mvnrnd(mup,Sigmap);


            Wy = y(pd+1:end);
            Sigmaq = inv(X'*diag(Wy)*X+inv(B));
            muq = Sigmaq*(X'*(Y-1/2)+inv(B)*b);
            betay = zeros(pd,1);

            qdpxx = normqdp(betax(:),mup,muq,Sigmap,Sigmaq);
            if rand <= qdpxx
                betay = betax(:);
            else
                while true
                    betayy = mvnrnd(muq,Sigmaq);
                    pdqyy = normqdp(betayy(:),muq,mup,Sigmaq,Sigmap);
                    if rand > pdqyy
                        betay = betayy(:);
                        break
                    end
                end
            end

            xk(1:pd) = betax;
            yk(1:pd) = betay';

            Zp = abs(X*betax(:));
            Wx = pgdraw(Zp);

            Zq = abs(X*betay(:));
            qdpxx = polygammaqdp(Zp,Zq,Wx);
            W = rand(n,1);
            Iac = W <= qdpxx;
            Wy = zeros(n,1);
            Wy(Iac) = Wx(Iac);

            Isampled = Iac;
            pdqyy = zeros(n,1);
            while(~all(Isampled))
                Wy(~Isampled) = pgdraw(Zq(~Isampled));
                pdqyy(~Isampled) = polygammaqdp(Zq(~Isampled),Zp(~Isampled),Wy(~Isampled));
                W(~Isampled) = rand(sum(~Isampled),1);
                Iac = Isampled | W > pdqyy;
                Isampled = Iac;
            end
            xk(pd+1:end) = Wx';
            yk(pd+1:end) = Wy';
        end
    end
end

function y = normqdp(x,mup,muq,sigmap,sigmaq)
invp = inv(sigmap);
invq = inv(sigmaq);
ap = mup'*invp*mup;
aq = muq'*invq*muq;
bq = 2*muq'*invq*x;
bp = 2*mup'*invp*x;
c = x'*(invp-invq)*x;
logqdp = log(det(sigmap))-log(det(sigmaq))+c+ap-bp+bq-aq;
y = exp(logqdp/2);
end

function y = polygammaqdp(Zp,Zq,Wx)
aq = log(cosh(Zq/2));
ap = log(cosh(Zp/2));
cqp = Zq.^2/2-Zp.^2/2;
y = exp(aq-ap-cqp.*Wx);
end