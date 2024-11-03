classdef model_probit
    %MODELCLASS Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ModelName      % Model name
        dimension
        dataY
        dataX
        XTX
        IXTXX
        IXTX
        SIXTX
        datalength
        init_sample
    end

    methods
        function obj = model_probit(dimension,Y,X,varargin)
            obj.ModelName  = 'probitGibbs';
            obj.dimension = dimension;
            obj.dataY = Y;
            obj.dataX = X;
            obj.datalength = length(X);
            obj.XTX = X'*X;
            obj.IXTX = inv(obj.XTX);
            obj.IXTXX = obj.IXTX*X';
            obj.SIXTX = chol(obj.IXTX,'lower');
            if nargin > 3
                obj.init_sample = sample0;
            else
                obj.init_sample = obj.initsample();
            end
        end

        %% Initialize parameters
        function x0 = initsample(obj)
            d = obj.dimension;
            X = obj.dataX;
            B = obj.IXTXX;
            Y = obj.dataY;
            beta0 = B*Y;
            z0 = normrnd(X*beta0,1);
            x0 = reshape([beta0;z0],1,d);
        end

        %% proposal sampler
        function y = proposal_sample(obj,x,u)
            d = obj.dimension;
            X = obj.dataX;
            C = obj.IXTXX;
            Y = obj.dataY;
            L = obj.SIXTX;
            n = obj.datalength;
            db = d-n;
            z = x(db+1:d);
            mu = C*z(:);
            u1 = reshape(u(1:db),[],1);
            u2 = reshape(u(db+1:end),[],1);
            beta = mu+L*norminv(u1);
            z = zeros(n,1);
            Y1 = Y == 1;
            Y0 = Y == 0;
            xb = X*beta(:);
            f0 = normcdf(-xb);
            z(Y1) = xb(Y1)+norminv(f0(Y1)+(1-f0(Y1)).*u2(Y1));
            z(Y0) = xb(Y0)+norminv(f0(Y0).*u2(Y0));
            y = reshape([beta;z],1,d);
        end

        %% proposal pdf
        function [xk,yk] = proposal_coupling_sample(obj,x,y,u)
            d = obj.dimension;
            X = obj.dataX;
            Isigma = obj.XTX;
            C = obj.IXTXX;
            Y = obj.dataY;
            L = obj.SIXTX;
            Sigma = obj.IXTX;
            n = obj.datalength;
            db = d-n;
            xk = zeros(1,d);
            yk = zeros(1,d);
            pz = x(db+1:d);
            qz = y(db+1:d);
            mup = C*pz(:);
            muq = C*qz(:);
            u1 = reshape(u(1:db),[],1);
            betax = mup+L*norminv(u1);%列向量
            betay = zeros(db,1);
            logqdpx = 1/2*(muq-mup)'*Isigma*(2*betax-muq-mup);
            if rand <= exp(logqdpx)
                betay = betax;
            else
                while true
                    byy = mvnrnd(muq,Sigma);
                    logpdqy = 1/2*(mup-muq)'*Isigma*(2*byy(:)-mup-muq);
                    if rand > exp(logpdqy)
                        betay = byy(:);
                        break
                    end
                end
            end
            xk(1:db) = betax';
            yk(1:db) = betay';
            Xbx = X*betax;
            Xby = X*betay;
            zx = zeros(n,1);
            qdpx = zeros(n,1);
            u2 = reshape(u(db+1:end),[],1);
            Y1 = Y == 1;
            Y0 = Y == 0;
            f0p = normcdf(-Xbx);
            f0q = normcdf(-Xby);
            muqmip = Xby-Xbx;
            muqplp = Xby+Xbx;
            zx(Y1) = Xbx(Y1)+norminv(f0p(Y1)+(1-f0p(Y1)).*u2(Y1));
            zx(Y0) = Xbx(Y0)+norminv(f0p(Y0).*u2(Y0));
            qdpx(Y1) = (1-f0p(Y1))./(1-f0q(Y1)).*exp(1/2*muqmip(Y1).*(2*zx(Y1)-muqplp(Y1)));
            qdpx(Y0) = f0p(Y0)./f0q(Y0).*exp(1/2*muqmip(Y0).*(2*zx(Y0)-muqplp(Y0)));
            W = rand(n,1);
            Iac = W <= qdpx;
            zy = zeros(n,1);
            zy(Iac) = zx(Iac);
            Isampled = Iac;
            pdqyy = zeros(n,1);
            while(~all(Isampled))
                u2(~Isampled) = rand(sum(~Isampled),1);
                IY1 = ~Isampled & Y == 1;
                IY0 = ~Isampled & Y == 0;
                zy(IY1) = Xby(IY1)+norminv(f0q(IY1)+(1-f0q(IY1)).*u2(IY1));
                zy(IY0) = Xby(IY0)+norminv(f0q(IY0).*u2(IY0));
                pdqyy(IY1) = (1-f0q(IY1))./(1-f0p(IY1)).*exp(1/2*(-muqmip(IY1)).*(2*zy(IY1)-muqplp(IY1)));
                pdqyy(IY0) = f0q(IY0)./f0p(IY0).*exp(1/2*(-muqmip(IY0)).*(2*zy(IY0)-muqplp(IY0)));
                W(~Isampled) = rand(sum(~Isampled),1);
                Iac = Isampled | W > pdqyy;
                Isampled = Iac;
            end
            xk(db+1:end) = zx';
            yk(db+1:end) = zy';
        end
    end

end

