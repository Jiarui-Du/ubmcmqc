classdef model_linear
    %MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        dataY
        dataX
        dataXTX
        dataXTY
        b0
        B0
        n0
        s0
        init_sample
    end

    methods
        function obj = model_linear(dimension,datay,dataX,b0,B0,n0,s0,varargin)
            obj.ModelName  = 'linearGibbs';
            obj.dimension = dimension;
            obj.dataY = datay;
            obj.dataX = dataX;
            obj.dataXTX = dataX'*dataX;
            obj.dataXTY = dataX'*datay;
            obj.b0 = b0;
            obj.B0 = B0;
            obj.n0 = n0;
            obj.s0 = s0;
            if nargin > 7
                obj.init_sample = sample0;
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function x0 = initsample(obj)
            d = obj.dimension;
            beta = mvnrnd(obj.b0,obj.B0);
            tau = 1/gamrnd(obj.n0/2,2/obj.s0);
            x0 = reshape([beta(:);tau],1,d);
        end

        %% proposal sampler
        function y = proposal_sample(obj,x,u)
            d = obj.dimension;
            Y = obj.dataY;
            X = obj.dataX;
            n = length(X);
            tau = x(end);
            invB0 = inv(obj.B0);
            invB0b = invB0*obj.b0;

            invtau = 1/tau;
            invB1 = invB0+invtau*obj.dataXTX;
            B1 = inv(invB1);
            b1 = B1*(invB0b+invtau*obj.dataXTY);
            L = chol(B1,'lower');
            u1 = reshape(u(1:d-1),[],1);

            beta = b1+L*norminv(u1);
            n1 = n + obj.n0;
            res = (Y-X*beta)'*(Y-X*beta);
            s1 = obj.s0+res;
            u2 = u(d);
            tau = 1/gaminv(u2,n1/2,2/s1);

            y = reshape([beta(:);tau],1,d);
        end

        %% proposal pdf
        function [xk,yk] = proposal_coupling_sample(obj,x,y,u)
            d = obj.dimension;
            X = obj.dataX;
            Y = obj.dataY;
            XTX = obj.dataXTX;
            XTY = obj.dataXTY;
            n = length(X);
            invBB = inv(obj.B0);
            bb = obj.b0;
            ss = obj.s0;
            nn = obj.n0;

            xk = zeros(1,d);
            yk = zeros(1,d);

            taup = x(end);
            itaup = 1/taup;
            B1p = inv(invBB+itaup*XTX);
            b1p = B1p*(invBB*bb+itaup*XTY);
            
            Lp = chol(B1p,'lower');
            u1 = reshape(u(1:d-1),[],1);
            betax = b1p+Lp*norminv(u1);

            tauq = y(end);
            itauq = 1/tauq;
            B1q = inv(invBB+itauq*XTX);
            b1q = B1q*(invBB*bb+itauq*XTY);

            qdpxx = normqdp(betax,b1p,b1q,B1p,B1q);
            if rand <= qdpxx
                betay = betax;
            else
                while true
                    byy = mvnrnd(b1q,B1q);
                    pdqyy = normqdp(byy(:),b1q,b1p,B1q,B1p);
                    if rand > pdqyy
                        betay = byy(:);
                        break
                    end
                end
            end

            xk(1:d-1) = betax';
            yk(1:d-1) = betay';

            n1 = n + nn;
            s1p = ss + (Y-X*betax)'*(Y-X*betax);
            u2 = u(d);
            taux = 1/gaminv(u2,n1/2,2/s1p);
            
            s1q = ss + (Y-X*betay)'*(Y-X*betay);
            qdpxx = invgammaqdp(taux,n1/2,s1p/2,s1q/2);

            if rand <= qdpxx
                tauy = taux;
            else
                while true
                    tauyy = 1/gamrnd(n1/2,2/s1q);
                    pdqyy = invgammaqdp(tauyy,n1/2,s1q/2,s1p/2);
                    if rand > pdqyy
                        tauy = tauyy;
                        break
                    end
                end
            end

            xk(d) = taux;
            yk(d) = tauy;
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

function y = invgammaqdp(x,alpha,bp,bq)
logqdp = alpha*(log(bq)-log(bp))-(bq-bp)/x;
y = exp(logqdp);
end
