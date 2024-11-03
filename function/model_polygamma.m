classdef model_polygamma
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
        samplemethod
    end

    methods
        function obj = model_polygamma(dimension,para_dimension,Y,X,b,B,method,varargin)
            obj.ModelName  = 'polygammaGibbs';
            obj.dimension = dimension;
            obj.para_dimension = para_dimension;
            obj.dataY = Y;
            obj.dataX = X;
            obj.dataB = B;
            obj.datab = b;
            obj.samplemethod = method;
            if nargin > 7
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
        function y = proposal_sample(obj,x,u)
            method = obj.samplemethod;
            d = obj.dimension;
            pd = obj.para_dimension;
            Y = obj.dataY;
            X = obj.dataX;
            B = obj.dataB;
            b = obj.datab;

            W = x(pd+1:end);
            Isigma = X'*diag(W)*X+inv(B);
            Sigma = inv(Isigma);
            mu = Sigma*(X'*(Y-1/2)+inv(B)*b);
            L = chol(Sigma,"lower");
            u1 = u(1:pd);
            beta = mu+L*norminv(u1(:));

            Z = abs(X*beta);
            
            if method == "mix1"
                u2 = u(pd+1:end);
                W = pgdraw_cudiidz(u2(:),Z);
            elseif method == "mix2"
                u2 = u(pd+1:2:end-1); %judge
                u3 = u(pd+2:2:end); % sample
                W = pgdraw_cudiidzu(u2(:),u3(:),Z);
            else
                u2 = u(pd+1:end);
                W = pgdrawnew(u2(:),Z);
            end
            y = reshape([beta(:);W(:)],1,d);
        end

        %% proposal pdf
        function [xk,yk] = proposal_coupling_sample(obj,x,y,u)
            method = obj.samplemethod;

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
            Lp = chol(Sigmap,"lower");
            u1 = u(1:pd);
            betax = mup+Lp*norminv(u1(:));


            Wy = y(pd+1:end);
            Sigmaq = inv(X'*diag(Wy)*X+inv(B));
            muq = Sigmaq*(X'*(Y-1/2)+inv(B)*b);
            Lq = chol(Sigmaq,"lower");
            betay = zeros(pd,1);

            qdpxx = normqdp(betax,mup,muq,Sigmap,Sigmaq);
            if rand <= qdpxx
                betay = betax;
            else
                while true
                    betayy = muq+Lq*norminv(rand(pd,1));
                    pdqyy = normqdp(betayy,muq,mup,Sigmaq,Sigmap);
                    wy = rand;
                    if wy > pdqyy
                        betay = betayy;
                        break
                    end
                end
            end

            xk(1:pd) = betax';
            yk(1:pd) = betay';

            Zp = abs(X*betax);
            if method == "mix1"
                u2 = reshape(u(pd+1:end),[],1);
                Wx = pgdraw_cudiidz(u2,Zp);
            elseif method == "mix2"
                u2 = reshape(u(pd+1:2:end-1),[],1); %judge
                u3 = reshape(u(pd+2:2:end),[],1); % sample
                Wx = pgdraw_cudiidzu(u2,u3,Zp);
            else
                u2 = reshape(u(pd+1:end),[],1);
                Wx = pgdrawnew(u2,Zp);
            end


            Zq = abs(X*betay);
            qdpxx = polygammaqdp(Zp,Zq,Wx);
            W = rand(n,1);
            Iac = W <= qdpxx;
            Wy = zeros(n,1);
            Wy(Iac) = Wx(Iac);

            Isampled = Iac;
            pdqyy = zeros(n,1);

            while(~all(Isampled))
                if method == "mix1"
                    u2(~Isampled) = rand(sum(~Isampled),1);
                    Wy(~Isampled) = pgdraw_cudiidz(u2(~Isampled),Zq(~Isampled));
                elseif method == "mix2"
                    u2(~Isampled) = rand(sum(~Isampled),1);
                    u3(~Isampled) = rand(sum(~Isampled),1);
                    Wy(~Isampled) = pgdraw_cudiidzu(u2(~Isampled),u3(~Isampled),Zq(~Isampled));
                else
                    u2(~Isampled) = rand(sum(~Isampled),1);
                    Wy(~Isampled) = pgdrawnew(u2(~Isampled),Zq(~Isampled));
                end
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