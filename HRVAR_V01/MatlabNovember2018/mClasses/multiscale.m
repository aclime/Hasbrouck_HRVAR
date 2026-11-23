classdef multiscale < matlab.mixin.Copyable
    %multiscale builds and contains synthesized multiscale irfs
    properties
        secondFactors
        irfms
        irfmst
    end
    methods
        function ms=multiscale(mva)
            secondFactors
            if nargin==0; return; end
            ms.secondFactor = mv.secondFactor;
            ms.irfms = mv.irf;
            ms.irfmst = mv.irft;
        end
        function build(ms, ms0, mv)
            %   ms0 has finer resolution than ms.
            ms.secondFactor = mv.secondFactor;
            ms.irfms = zeros(size(mv.irf));
            ms.irfmst = mv.irft;
            k = find(ms0.irfmst==10,1);
            u = squeeze( ms0.irfms(:,k,:) );
            for k=1:size(mv.irf,2)
                ms.irfms(:,k,:) = squeeze(mv.irf(:,k,:))*u;
            end
            ms.irfmst = mv.irft+1;
        end
        function buildAlt(ms,ms0,mv)    %   not used!!!!!
            %  Direct IRF computation
            nAhead = mv.irft(end);
            ms.secondFactor = mv.secondFactor;
            if mv.arP>0
                phi=flip(mv.phiSet,3);
                phi=reshape(phi,mv.nPrices,mv.nPrices*mv.arP);
            end
            irf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            dp1 = squeeze(ms0.irfms(:,11,:));
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,mv.arP+nAhead);
                dp(:,1+mv.arP) = dp1(:,iShock);
                p = dp1(:,iShock);
                for t=(2+mv.arP):(nAhead+mv.arP)
                    d = zeros(mv.nPrices,1);
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = d + gamma*B*p; end;
                    dp(:,t) = d;
                    if mv.arP>0; dp(:,t) = d + phi*reshape(dp(:,t-mv.arP:t-1),mv.nPrices*mv.arP,1); end
                    p = p + dp(:,t);
                end
%                 if display
%                     nPrint = min(nAhead,10);
%                     fprintf('MVARf.irfBuild iShock=%d\n',iShock)
%                     disptable(dp(:,1:10),strsplit(int2str(1:10)),mv.priceNames);
%                 end;
                irf(:,:,iShock) = cumsum(dp(:,mv.arP:end),2);
            end
            ms.irfms = irf;
            ms.irfms(:,1:11,:);
            ms.irfmst = 0:size(irf,2)-1;
%             if display
%                 for i=1:mv.nPrices;
%                     x = [irf(:,1:nPrint,i) irf(:,end,i)];
%                     s = strsplit(int2str([0:nPrint-1 size(irf,2)-1]));
%                     disptable(x,s,mv.priceNames);
%                 end
%             end
        end
    end
end

