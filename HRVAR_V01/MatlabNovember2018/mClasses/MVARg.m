classdef MVARg < matlab.mixin.Copyable
    % Microstructure VAR/VECM class
    % Contains for both autoregressive (AR) and PDL specifications.
    % Only the PDL are checked.
    % The notional data matrix is z, and the notional crossproduct matrix is z'p (zpz).
    % The column names of z are set in setnamesPDL or setnamesAR
    % [1 (P1 block) (P2 block) ... (PnPrice block) eVecs]
    % where a price block includes the price (not lagged) and all polynomials or AR coefficients.
    % the eVecs are the error correction terms.
    properties
        symbol  %   ticker symbol
        inDate   %   sample date
        showTiming = true
        secondFactor = 1        % The number of time units per second.
        prices=spl.empty
        dPrices=spd.empty
        nPrices
        %   Specification variables
        spec=Specification.PDL;         %   specification type: AR or PDL
        intercept=true  %   include intercept in specifications
        ecm=true        %   include error correction terms
        firstValid %   first valid time (index) in the 'data' matrix
        lastValid  %   last ...
        %   if spec=AR, then arP must be set; if spec=PDL then polys must be set.
        arP;                    % order of VAR
        polys=polynom.empty     % array of polynomials that define the design matrix for the coefficients
        kzlhs   %   logical index vector within z to left hand side vars
        kzrhs   %   ... and right hand side vars
        priceNames
        zNames
        xNames
        yNames
        eVecs=spl.empty   %   error (cointegration) vectors
        eVecMeans
        zpz    %   full cross-product matrix
        %   parameters
        b      %   coefficients
        bRowNames
        bColNames
        seb     %   se's of coefficients
        tb      %   t's of coefficients
        eCov    %   residual covariance matrix
        eCorr   %   residual correlation matrix
        irft    %   time indices for irf array (intervals; to get to real time, divide by secondFactor)
        irf     %   nPrices x nAhead x nPrices(=nShocks)
        irfPacked   %   Condensed
        bIRFt   %   time indices for bridged irf (intervals; to get to real time, divide by secondFactor)
        bIRF    %   bridged IRF 
        rwd     %   random walk decomposition (object)
    end
    methods
        %% initialization and setup routines.
        function o=MVARg(n) % array constructor
            %   MVARg(n) returns an array of MVARg's.
            if nargin~=0; o(n)=MVARg; end
        end
        function setup(mv,mDisplay)
            %   take price differences, form cointegrating vectors, etc.
            if nargin<2; mDisplay=false; end
            if mDisplay; fprintf('setup.\n'); end
            mv.nPrices = length(mv.prices);
            mv.dPrices(mv.nPrices) = spd;
            for i=1:mv.nPrices
                if mv.prices(i).maxSize~=mv.prices(1).maxSize; error('MVARg.setup. price vectors don''t have the same maxSize.'); end
                mv.dPrices(i) = mv.prices(i).diff;
            end
            if mv.ecm && mv.nPrices>=2
                if mDisplay; fprintf('Computing eVecs\n'); end
                mv.eVecs(mv.nPrices-1)=spl;
                for i=1:mv.nPrices-1
                    %fprintf('%s [minus]\n',mv.prices(1).toString)
                    %fprintf('%s [equals]\n',mv.prices(i+1).toString)
                    mv.eVecs(i) = mv.prices(1).splMinus(mv.prices(i+1));
                    mv.eVecs(i).name = [mv.prices(1).name '-' mv.prices(i+1).name];
                    if mDisplay; fprintf('%s\n',mv.eVecs(i).toString); end
                end
            end
        end
        function setNamesAR(mv,nDisplay)   %   AR specification set up row and column labels for various matrices;
            %   set column names of (implicit) data matrix. nDisplay=true to display the name arrays.
            if nargin==1; nDisplay=false; end
            mv.priceNames = cellstr( arrayfun(@(x) x.name,mv.prices,'UniformOutput',false) );
            mv.zNames = {'const'};
            for i=1:mv.nPrices
                mv.zNames=horzcat(mv.zNames, mv.dPrices(i).name);
                if mv.arP>0; mv.zNames = horzcat(mv.zNames, arrayfun(@(x) strcat(mv.dPrices(i).name, "L", int2str(x)), 1:mv.arP)); end
            end
            if mv.ecm
                for i=1:mv.nPrices-1
                    mv.zNames = [mv.zNames mv.eVecs(i).name];
                end
            end
            k = [false repmat([true, false(1,mv.arP)],1,mv.nPrices)];
            if mv.ecm; k=[k false(1,mv.nPrices-1)]; end
            mv.kzlhs = k;
            mv.yNames = mv.zNames(mv.kzlhs);
            mv.kzrhs = [mv.intercept repmat([false true(1,mv.arP)],1,mv.nPrices)];
            if mv.ecm; mv.kzrhs = [mv.kzrhs true(1,mv.nPrices-1)]; end
            mv.xNames = mv.zNames(mv.kzrhs);
            mv.zNames = cellstr(mv.zNames);
            mv.yNames = cellstr(mv.yNames);
            mv.xNames = cellstr(mv.xNames);
            if nDisplay
                disp('zNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.zNames))
                if mod(length(mv.zNames),6)~=0; fprintf('\n'); end
                disp('yNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.yNames))
                if mod(length(mv.yNames),6)~=0; fprintf('\n'); end
                disp('xNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.xNames))
                if mod(length(mv.xNames),6)~=0; fprintf('\n'); end
                %disp(['kzlhs: ' int2str(mv.kzlhs)])
                %disp(['kzrhs: ' int2str(mv.kzrhs)])
            end
            %mv.kpSelect = true(1,mv.polys.nTerms);  %   default; may be overridden.
        end
        function setNamesPDL(mv,nDisplay)   %   set up row and column labels for various matrices;
            %   set column names of (implicit) data matrix. nDisplay=true to display the name arrays.
            if nargin==1; nDisplay=false; end
            pn = {''};
            for i=1:length(mv.polys); pn = horzcat(pn, mv.polys(i).vNames); end
            mv.priceNames = cellstr( arrayfun(@(x) x.name,mv.prices,'UniformOutput',false) );
            mv.zNames = {'const'};
            for i=1:mv.nPrices
                mv.zNames = horzcat(mv.zNames, strcat(mv.dPrices(i).name,pn));
            end
            if mv.ecm
                for i=1:mv.nPrices-1
                    mv.zNames = [mv.zNames mv.eVecs(i).name];
                end
            end
            m = length(pn)-1;
            k = [false repmat([true, false(1,m)],1,mv.nPrices)];
            if mv.ecm; k=[k false(1,mv.nPrices-1)]; end
            mv.kzlhs = k;
            mv.yNames = mv.zNames(mv.kzlhs);
            mv.kzrhs = [mv.intercept repmat([false true(1,m)],1,mv.nPrices)];
            if mv.ecm; mv.kzrhs = [mv.kzrhs true(1,mv.nPrices-1)]; end
            mv.xNames = mv.zNames(mv.kzrhs);
            mv.zNames = cellstr(mv.zNames);
            mv.yNames = cellstr(mv.yNames);
            mv.xNames = cellstr(mv.xNames);
            if nDisplay
                disp('zNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.zNames))
                if mod(length(mv.zNames),6)~=0; fprintf('\n'); end
                disp('yNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.yNames))
                if mod(length(mv.yNames),6)~=0; fprintf('\n'); end
                disp('xNames:');
                fprintf([repmat('%-18s ',1,6) '\n'],string(mv.xNames))
                if mod(length(mv.xNames),6)~=0; fprintf('\n'); end
                %disp(['kzlhs: ' int2str(mv.kzlhs)])
                %disp(['kzrhs: ' int2str(mv.kzrhs)])
            end
        end
        function setValid(mv,maxLag,firstTime,lastTime,mDisplay)
            % setValid(mv,maxLag,firstTime,lastTime,mDisplay)
            % firstTime and lastTime are duration(h,m,s)
            if nargin<2; mDisplay=false; end
            if mDisplay; fprintf('setValid (computing largest possible range)\n'); end
            %first=0;
            maxSize = mv.prices(1).maxSize;
            first=seconds(firstTime)*mv.secondFactor;
            spd.setgetMax(maxSize);
            %last = maxSize+1;
            last = seconds(lastTime)*mv.secondFactor;
            timePrec = log10(mv.secondFactor);
            %last=spd.setgetMax+1;
            %             maxLag = mv.polys.maxLength;
            for i=1:mv.nPrices
                tFirst = timeSec(mv.prices(i).firstValid/mv.secondFactor,timePrec);
                tLast =  timeSec(mv.prices(i).lastValid/mv.secondFactor,timePrec);
                if mDisplay; fprintf('%12s %12d (%s) to %12d (%s)\n',mv.prices(i).name, mv.prices(i).firstValid, tFirst, mv.prices(i).lastValid, tLast); end
                if mv.prices(i).maxSize~=maxSize; error('price vectors don''t have the same maxSize.'); end
                first=max([first,mv.prices(i).firstValid+maxLag+1]);
                last=min([last,mv.prices(i).lastValid]);
                % if mDisplay; fprintf('first=%d last=%d\n',first,last); end
            end
            tFirst = timeSec(first/mv.secondFactor,timePrec);
            tLast =  timeSec(last/mv.secondFactor,timePrec);
            mv.firstValid = first;
            mv.lastValid = last;
            fprintf('valid range: %12d (%s) to %12d (%s)\n',mv.firstValid, tFirst, mv.lastValid, tLast);
        end
        function eVecDemean(mv,eDisplay)
            if ~mv.ecm; return; end
            if nargin<2; eDisplay=false; end
            mv.eVecMeans = NaN(1,mv.nPrices-1);
            n = mv.lastValid - mv.firstValid + 1;
            for i=1:mv.nPrices-1
                z = mv.eVecs(i).sum(mv.firstValid, mv.lastValid, 1);
                mv.eVecMeans(i) = z/n;
                mv.eVecs(i).v = mv.eVecs(i).v - mv.eVecMeans(i);
            end
            if eDisplay; disp(['eVecMeans: ' num2str(mv.eVecMeans)]); end
        end
        function pricePlot(mv)
            %   plot from (sparse) input prices
            mv.nPrices = length(mv.prices);
            T = spd.setgetMax();
            maxlen = 0;
            for i=1:mv.nPrices
                maxlen = max([maxlen,length(mv.prices(i).i)]);
            end
            X = T*ones(maxlen,mv.nPrices);
            Y = zeros(maxlen,mv.nPrices);
            for i=1:mv.nPrices
                pLast = mv.prices(i).v(end);
                Y(:,i) = pLast;
                m = length(mv.prices(i).i);
                X(1:m,i)=mv.prices(i).i;
                Y(1:m,i)=mv.prices(i).v;
            end
            stairs(X,Y);
        end
        %% IRF routines
        function irf=irfBuild2R(mv,nAhead,iDisplay) %  Direct IRF computation
            if nargin<3; iDisplay=false; end
            phi=mv.phiSet;
            maxLag = size(phi,3);    % This may be <= the maxLag in the zpz matrix. ????? not in this version!
            irf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            T = maxLag + nAhead+1;
            dpf = zeros(mv.nPrices,nAhead+1);
            mPhi = reshape(phi,mv.nPrices,[]);   %   Restack phi into an nPrice x (nPrice*maxLag) maxtrix
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,T);
                p = zeros(mv.nPrices,1);
                e = zeros(mv.nPrices,1);
                e(iShock) = 1;
                for t=T-maxLag:-1:1
                    %fprintf('i=%d t=%d\n',i,t)
                    d = 0;
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = gamma*B*p; end;
                    v = reshape(dp(:,t+1:t+maxLag),[],1);
                    %disp(['size(v)=' int2str(size(v))])
                    dp(:,t) = d+mPhi*v+e;
                    p = p + dp(:,t);
                    e = zeros(mv.nPrices,1);
                end
                dpf(:,:,iShock) = flip(dp(:,1:nAhead+1),2);
                irf(:,:,iShock) = cumsum( flip(dp(:,1:nAhead+1),2), 2);
            end
            mv.irf = irf;
            if iDisplay
%                 disp('irfBuild2R dpIRF:')
%                 for i=1:mv.nPrices; disptableLog(dpf(:,:,i),mv.priceNames,5,0); end
                disp('irfBuild2R pIRF:')
                for i=1:mv.nPrices; disptableLog(irf(:,:,i),mv.priceNames,5,0); end
            end
        end
        function irf=irfBuild2(mv,nAhead,iDisplay) %  (UPDATED 7/18/18) Direct IRF computation
            if nargin<3; iDisplay=false; end
            phi = mv.phiSet;
            maxLag = size(phi,3);
            T = maxLag + nAhead + 1;
            irf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            dpf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            mPhi = reshape(flip(phi,3), mv.nPrices,[]);
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,T);
                p = zeros(mv.nPrices,1);
                e0 = zeros(mv.nPrices,1);
                e0(iShock) = 1;
                for t=maxLag+1:T
                    d = 0;
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = gamma*B*p; end;
                    v = reshape(dp(:,t-maxLag:t-1),[],1);
                    %v = reshape(dp(:,t+1:t+maxLag),[],1);
                    %disp(['size(v)=' int2str(size(v))])
                    dp(:,t) = d+mPhi*v;
                    if t==maxLag+1; dp(:,t) = dp(:,t)+e0; end
                    p = p + dp(:,t);
                end
                dpf(:,:,iShock) = dp(:,maxLag+1:T);
                irf(:,:,iShock) = cumsum(dp(:,maxLag+1:T),2);
            end
            mv.irf = irf;
            if iDisplay
%                 disp('irfBuild2 dpIRF:')
%                 for i=1:mv.nPrices; disptableLog(dpf(:,:,i),mv.priceNames,5,0); end
                disp('irfBuild2 pIRF:')
                for i=1:mv.nPrices; disptableLog(irf(:,:,i),mv.priceNames,5,0); end
            end
        end
        function [irfPacked,irft]=irfBuild2Packed(mv,nAhead,iDisplay) %  (UPDATED 7/18/18) Direct IRF computation
            if nargin<3; iDisplay=false; end
            phi = mv.phiSet;
            maxLag = size(phi,3);
            T = maxLag + nAhead + 1;
            m = ceil(log10(nAhead));
            k = unique( horzcat(reshape( (0:10)'*(10.^(0:m)), 1, []), nAhead) );
            irft = k(k<=nAhead);
            % s = strsplit(int2str(k));
            irf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            dpf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
            mPhi = reshape(flip(phi,3), mv.nPrices,[]);
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,T);
                p = zeros(mv.nPrices,1);
                e0 = zeros(mv.nPrices,1);
                e0(iShock) = 1;
                for t=maxLag+1:T
                    d = 0;
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = gamma*B*p; end;
                    
                    v = reshape(dp(:,t-maxLag:t-1),[],1);
                    %v = reshape(dp(:,t+1:t+maxLag),[],1);
                    %disp(['size(v)=' int2str(size(v))])
                    dp(:,t) = d+mPhi*v;
                    if t==maxLag+1; dp(:,t) = dp(:,t)+e0; end
                    p = p + dp(:,t);
                end
                dpf(:,:,iShock) = dp(:,maxLag+1:T);
                irf(:,:,iShock) = cumsum(dp(:,maxLag+1:T),2);
            end
            dpfPacked = dpf(:,irft+1,:);
            irfPacked = irf(:,irft+1,:);
            if iDisplay
%                 disp('irfBuild2Packed dpIRF:')
%                 for i=1:mv.nPrices; disptable(dpfPacked(:,:,i),strsplit(int2str(irft)),mv.priceNames); end
                disp('irfBuild2Packed pIRF:')
                for i=1:mv.nPrices; disptable(irfPacked(:,:,i),strsplit(int2str(irft)),mv.priceNames); end
            end
        end
        function irfBuildAR(mv,nAhead,rDisplay) %  Direct IRF computation
            if nargin<3; rDisplay=false; end
            if mv.arP>0
                phi=flip(mv.phiSetAR,3);
                phi=reshape(phi,mv.nPrices,mv.nPrices*mv.arP);
            end
            irf = zeros(mv.nPrices,nAhead,mv.nPrices);
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,mv.arP+nAhead);
                e = zeros(mv.nPrices,1);
                e(iShock) = 1;
                p = e;
                dp(:,1+mv.arP) = e;
                for t=(2+mv.arP):(nAhead+mv.arP)
                    d = zeros(mv.nPrices,1);
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = d + gamma*B*p; end;
                    dp(:,t) = d;
                    if mv.arP>0; dp(:,t) = d + phi*reshape(dp(:,t-mv.arP:t-1),mv.nPrices*mv.arP,1); end
                    p = p + dp(:,t);
                end
                if rDisplay
                    nPrint = min(nAhead,10);
                    fprintf('MVARg.irfBuild iShock=%d\n',iShock)
                    disptable(dp(:,1:10),strsplit(int2str(1:10)),mv.priceNames);
                end;
                irf(:,:,iShock) = cumsum(dp(:,1+mv.arP:end),2);
            end
            mv.irf = irf;
            mv.irft = 0:size(irf,2)-1;
            if rDisplay
                for i=1:mv.nPrices;
                    x = [irf(:,1:nPrint,i) irf(:,end,i)];
                    s = strsplit(int2str([0:nPrint-1 size(irf,2)-1]));
                    disptable(x,s,mv.priceNames);
                end
            end
        end
        function irfBuildPDL(mv,nAhead,nPrint) %  Direct IRF computation
            if nargin<3; nPrint=0; end
            phi=flip(mv.phiSetPDL,3);
            phi=reshape(phi,mv.nPrices,[]);
            maxLag = mv.polys.maxLength;
            irf = zeros(mv.nPrices,nAhead,mv.nPrices);
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
                B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            end
            for iShock=1:mv.nPrices  %   index of initial shock
                dp = zeros(mv.nPrices,maxLag+nAhead);
                e = zeros(mv.nPrices,1);
                e(iShock) = 1;
                p = e;
                dp(:,1+maxLag) = e;
                for t=(2+maxLag):(nAhead+maxLag)
                    d = zeros(mv.nPrices,1);
                    if mv.intercept; d = d + mv.b(1,:)'; end
                    if mv.ecm; d = d + gamma*B*p; end;
                    dp(:,t) = d + phi*reshape(dp(:,t-maxLag:t-1),mv.nPrices*maxLag,1);
                    p = p + dp(:,t);
                end
                irf(:,:,iShock) = cumsum(dp(:,1+maxLag:end),2);
            end
            mv.irf = irf;
            mv.irft = 0:size(irf,2)-1;
            if nPrint>0
                for i=1:mv.nPrices
                    disptableLog(irf(:,:,i),mv.priceNames,5,0);
                end
            end
        end
        function irf=irfBuildf(mv,nAhead,phi0,iDisplay) %   Direct IRF Computation with parfor construction
            %   There are nAhead forecasts; p0 is prepended to the beginning
            if nargin<4; iDisplay=false; end
            maxLag=size(phi0,3);
            irf = zeros(mv.nPrices, nAhead+1, mv.nPrices);
            parfor i=1:mv.nPrices
                pLag = zeros(mv.nPrices,1)
                pLag(i) = 1;
                f = mv.vForecast(pLag,phi0,nAhead);
                irf(:,:,i) = cumsum([pLag(:,end) f],2);
                %             f = p0+cumsum(dp(:,maxLag+1:end),2);
            end
            if iDisplay
                disp('irfBuildf irf:')
                for i=1:mv.nPrices
                    disptableLog(irf(:,:,i),mv.priceNames,5,0)
                end
            end
        end
        function f=vForecast(mv,pLag,phi0,nAhead,fDisplay)
            % pLag is ... p(t-2) p(t-1)
            % f is dp(t) ... dp(t+nAhead-1)
            if nargin<6; fDisplay=false; end;
            n = mv.nPrices;
            pLag = [zeros(n,1) pLag];  % pad with zeros so that can take differences.
            p = pLag(:,end);
            dpLag = diff(pLag,1,2);
            phi = reshape(flip(phi0,3),n,[]);
            intercept = 0;
            if mv.intercept; intercept=mv.b(1,:)'; end
            gammaB=0;
            if mv.ecm
                gamma = mv.b(size(mv.b,1)-n+2:end,:)';
                B = [ones(n-1,1) -eye(n-1)];
                gammaB = gamma*B;
            end
            maxLag = size(phi0,3);
            dp = zeros(n,nAhead+maxLag);
            nLags = size(dpLag,2);
            dp(:,maxLag-nLags+1:maxLag) = dpLag;
            for t=maxLag+1:nAhead+maxLag
                dp(:,t) = intercept + gammaB*p + phi*reshape(dp(:,t-maxLag:t-1),mv.nPrices*maxLag,1);
                p = p + dp(:,t);
            end
            f = dp(:,maxLag+1:end);
            if fDisplay; disptableLog(f,mv.priceNames,5); end
        end
        function irf=irfPDLr(mv,nAhead,fDisplay)
            if nargin<3; fDisplay=false; end
            irf = NaN(mv.nPrices,nAhead+1,mv.nPrices);
            parfor iShock=1:mv.nPrices
                e0 = zeros(mv.nPrices,1);
                e0(iShock)=1;
                f = fPDLr(mv,e0,nAhead+1,fDisplay);
                irf(:,:,iShock) = f;
            end
            if fDisplay
                disp('irfPDLr');
                for i=1:mv.nPrices
                    disptableLog(irf(:,:,i),mv.priceNames,5,0)
                end
            end
        end
        function f=fPDLr(mv,e0,nAhead,fDisplay)   %   forecast PDL, with state variables.
            if nargin<4; fDisplay=false; end
            maxLag = mv.polys.maxLength;
            T = maxLag + nAhead;
            gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
            B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            %   Extract coefficients of state vectors.
            k0 = 1;
            if mv.intercept; k0=2; end
            k1 = size(mv.b,1);
            if mv.ecm; k1 = k1 - mv.nPrices + 1; end
            bsc = mv.b(k0:k1,:);    % bsc: b state coefficients
            nPolys = length(mv.polys);
            dp = zeros(mv.nPrices,T);
            s = zeros(mv.polys.nTerms,mv.nPrices);    %   state vector
            p = zeros(mv.nPrices,1);
            for t=T-maxLag:-1:1
                d = 0;
                if mv.intercept; d = d + mv.b(1,:)'; end
                if mv.ecm
                    Bp = repmat(p(1),mv.nPrices-1,1) - p(2:end);
                    d = d + gamma*Bp;
                end
                %   update prices
                dp(:,t) = d + bsc' * reshape(s,[],1) + e0;
                p = p + dp(:,t);    %   cumulate price changes
                e0 = zeros(mv.nPrices,1);    %   except for the first impulse, all es are zero.
                %   update state vector
                s0 = s;
                yIn = zeros(nPolys,mv.nPrices);
                yOut = zeros(nPolys,mv.nPrices);
                for i=1:nPolys
                    yIn(i,:) = dp(:,t+mv.polys(i).kOffset);
                    yOut(i,:) = dp(:,t+mv.polys(i).kOffset+mv.polys(i).n);
                end
                s = mv.polys.stateUpdate(s0,yIn,yOut);
            end
            f = cumsum( flip(dp(:,1:nAhead),2), 2);
        end
        function irf=irfPDL(mv,nAhead,fDisplay)     % (UPDATED 7/18/18) 
            if nargin<3; fDisplay=false; end
            irf = NaN(mv.nPrices,nAhead+1,mv.nPrices);
            parfor iShock=1:mv.nPrices
                e0 = zeros(mv.nPrices,1);
                e0(iShock)=1;
                f = fPDL(mv,e0,nAhead,fDisplay);
                irf(:,:,iShock) = f;
            end
            if fDisplay
                disp('irfPDL');
                for i=1:mv.nPrices
                    disptableLog(irf(:,:,i),mv.priceNames,5,0)
                end
            end
        end
        function f=fPDL(mv,e0,nAhead,fDisplay)   %   (UPDATED 7/18/18) forecast PDL, with state variables.
            if nargin<4; fDisplay=false; end
            maxLag = mv.polys.maxLength;
            T = maxLag + nAhead + 1;
            gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
            B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            %   Extract coefficients of state vectors.
            k0 = 1;
            if mv.intercept; k0=2; end
            k1 = size(mv.b,1);
            if mv.ecm; k1 = k1 - mv.nPrices + 1; end
            bsc = mv.b(k0:k1,:);    % bsc: b state coefficients
            nPolys = length(mv.polys);
            dp = zeros(mv.nPrices,T);
            s = zeros(mv.polys.nTerms,mv.nPrices);    %   state vector
            p = zeros(mv.nPrices,1);
            for t=maxLag+1:T
                d = 0;
                if mv.intercept; d = d + mv.b(1,:)'; end
                if mv.ecm
                    Bp = repmat(p(1),mv.nPrices-1,1) - p(2:end);
                    d = d + gamma*Bp;
                end
                %   update prices
                dp(:,t) = d + bsc' * reshape(s,[],1);
                if t==maxLag+1; dp(:,t) = dp(:,t) + e0; end
                p = p + dp(:,t);    %   cumulate price changes
                %   update state vector
                s0 = s;
                yIn = zeros(nPolys,mv.nPrices);
                yOut = zeros(nPolys,mv.nPrices);
                for i=1:nPolys
                    yIn(i,:) = dp(:,t-mv.polys(i).kOffset);
                    yOut(i,:) = dp(:,t-(mv.polys(i).kOffset+mv.polys(i).n));
                end
                s = mv.polys.stateUpdate(s0,yIn,yOut);
            end
            f = cumsum(dp(:,1+maxLag:T),2);
          end
        function [irfPacked,irft]=irfPDLpacked(mv,nAhead,fDisplay)     % (UPDATED 7/18/18) 
            if nargin<3; fDisplay=false; end
            m = ceil(log10(nAhead));
            k = unique( horzcat(reshape( (0:10)'*(10.^(0:m)), 1, []), nAhead) );
            irft = k(k<=nAhead);
%             if fDisplay
%                 fprintf('irfPDLpacked. irft:\n')
%                 disptable(irft,split(int2str(1:length(irft))))
%             end
            irfPacked = NaN(mv.nPrices,length(irft),mv.nPrices);
            parfor iShock=1:mv.nPrices
%             for iShock=1:1
                e0 = zeros(mv.nPrices,1);
                e0(iShock)=1;
                f = fPDLpacked(mv,e0,nAhead,irft,fDisplay);
                irfPacked(:,:,iShock) = f;
            end
            if fDisplay
                disp('irfPDLpacked irfPacked');
                for i=1:mv.nPrices
                    disptable(irfPacked(:,:,i),split(int2str(irft)),mv.priceNames);
                end
            end
        end
        function fPacked=fPDLpacked(mv,e0,nAhead,irft,fDisplay)   %   (UPDATED 7/18/18) forecast PDL, with state variables, packed, buffering.
            %   irft is the row vector of time indices to save: 0, 1, 2, ..., nAhead
            if nargin<4; fDisplay=false; end
            maxLag = mv.polys.maxLength;
            T = maxLag + nAhead + 1;
            nShift = 10000;
            %   Tb defines a buffer for the forecasts. When the buffer is full, we shift, effectively discarding earlier values that aren't needed.
            Tb = min(T, maxLag+1+nShift);
            offset = 0;
            pShift = zeros(mv.nPrices,1);   %   price level immediately before the start of the buffer.
            dp = zeros(mv.nPrices,Tb);
            fPacked = zeros(mv.nPrices,size(irft,2));
            gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
            B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            %   Extract coefficients of state vectors.
            k0 = 1;
            if mv.intercept; k0=2; end
            k1 = size(mv.b,1);
            if mv.ecm; k1 = k1 - mv.nPrices + 1; end
            bsc = mv.b(k0:k1,:);    % bsc: b state coefficients
            nPolys = length(mv.polys);
            dp = zeros(mv.nPrices,T);
            s = zeros(mv.polys.nTerms,mv.nPrices);    %   state vector
            p = zeros(mv.nPrices,1);
            for t=maxLag+1:T
                td = t-offset;
                d = 0;
                if mv.intercept; d = d + mv.b(1,:)'; end
                if mv.ecm
                    Bp = repmat(p(1),mv.nPrices-1,1) - p(2:end);
                    d = d + gamma*Bp;
                end
                %   update prices
                dp(:,td) = d + bsc' * reshape(s,[],1);
                if t==maxLag+1; dp(:,td) = dp(:,td) + e0; end
                p = p + dp(:,td);    %   cumulate price changes
                %   update state vector
                s0 = s;
                yIn = zeros(nPolys,mv.nPrices);
                yOut = zeros(nPolys,mv.nPrices);
                for i=1:nPolys
                    yIn(i,:) = dp(:,td-mv.polys(i).kOffset);
                    yOut(i,:) = dp(:,td-(mv.polys(i).kOffset+mv.polys(i).n));
                end
                s = mv.polys.stateUpdate(s0,yIn,yOut);
                if t==T | td==Tb %  Save and shift.
%                     fprintf('Offset=%d\n',offset);
%                     disp('pShift:')
%                     disptable(pShift')
                    dpCum = pShift + cumsum(dp,2);
                    dpt = offset - maxLag + (0:Tb); %   time indices of values in buffer
                    [lia,locb] = ismember(irft,dpt);
                    fPacked(:,lia) = dpCum(:,locb(lia));
                    % disptable(fPacked,strsplit(int2str(irft)),mv.priceNames)
                    if t<T  % discard oldest nShift values
                        pShift = dpCum(:,nShift);
                        dp = circshift(dp,-nShift,2);
                        dp(:,maxLag+2:end) = NaN;
                        offset = offset + nShift;
                    end
                end
            end
          end
        %% constructing the zpz matrix
        function buildzpzSparseAR(mv) %   January, 2018
            %   build zpz using calls to sparse matrix routines.
            first = mv.firstValid;
            last = mv.lastValid;
            m = length(mv.zNames);
            zpz = zeros(m);
            zpz(1,1) = mv.lastValid-mv.firstValid+1;
            %   first column ...
            for i=1:mv.nPrices
                k = 1 + (i-1)*(1+mv.arP) + (1:mv.arP+1);
                zpzPart = Crossproduct.lagmxconst(mv.dPrices(i),mv.arP,first,last);
                zpz(k,1) = zpzPart;
            end
            %   blocks across prices
            for ip=1:mv.nPrices
                kr = 1 + (ip-1)*(mv.arP+1) + (1:mv.arP+1);
                for jp=1:ip
                    kc = 1 + (jp-1)*(mv.arP+1) + (1:mv.arP+1);
                    zpz(kr,kc) = Crossproduct.lagmxlagm(mv.dPrices(ip), mv.arP, mv.dPrices(jp), mv.arP, first, last);
                end
            end
            if mv.ecm
                for icv=1:length(mv.eVecs)
                    kr = 1 + (mv.arP + 1)*mv.nPrices + icv;
                    zpz(kr,1) = Crossproduct.splxconst(mv.eVecs(icv),1,first,last);
                    for jp=1:mv.nPrices
                        kc = 1 + (jp-1)*(mv.arP+1) + (1:mv.arP+1);
                        zpz(kr,kc) = Crossproduct.lagmxspl(mv.dPrices(jp), mv.arP, mv.eVecs(icv), 1, first, last);
                    end
                    for jcv=1:icv
                        kc = 1 + (mv.arP + 1)*mv.nPrices + jcv;
                        zpz(kr,kc) = Crossproduct.splxspl(mv.eVecs(icv), 1, mv.eVecs(jcv), 1, first, last);
                    end
                end
            end
            %   fill in upper triangle of zpz.
            zpz = tril(zpz,-1) + diag(diag(zpz)) + tril(zpz,-1)';
            mv.zpz = zpz;
        end
        function buildzpzSparsePDL(mv)
            %   build zpz using calls to sparse matrix routines.
            first = mv.firstValid;
            last = mv.lastValid;
            m = length(mv.zNames);
            zpz = zeros(m);
            zpz(1,1) = mv.lastValid-mv.firstValid+1;
            %   first column ...
            for i=1:mv.nPrices
                k = (i-1)*(1+mv.polys.nTerms) + 1;
                zpz(k+1,1) = Crossproduct.spdxconst(mv.dPrices(i),0,first,last);
                kk = k+2 : k+mv.polys.nTerms+1;
                zpz(kk,1) = Crossproduct.polyxconst(mv.dPrices(i),mv.polys,1,first,last);
            end
            %   blocks across prices
            nTerms = mv.polys.nTerms;
            dPrices = mv.dPrices.copy;
            polys = mv.polys.copy;
            nPrices = mv.nPrices;
            parfor iv=1:nPrices
                k = 1 + (iv-1)*(1+nTerms);
                kiv = (k+1):(k+nTerms+1);
                zpz0 = zeros(m);
                for jv=1:iv
                    jTimer = tic
                    fprintf('Starting iv=%d (%12s) jv=%d (%12s) ... ',iv,dPrices(iv).name,jv,dPrices(jv).name)
                    k = 1 + (jv-1)*(1+nTerms);
                    kjv = (k+1):(k+nTerms+1);
                    zpzij = zeros(1+nTerms);
                    zpzij(1,1) = Crossproduct.spdxspd(dPrices(iv), 0, dPrices(jv), 0, first, last);
                    zpzij(2:end,1) = Crossproduct.polyxspd(dPrices(iv), polys, 1, dPrices(jv), 0, first,last);
                    zpzij(1,2:end) = Crossproduct.polyxspd(dPrices(jv), polys, 1, dPrices(iv), 0, first, last);
                    zpzij(2:end,2:end) = Crossproduct.polyxpoly(dPrices(iv), polys, 1, dPrices(jv), polys, 1, first, last);
                    zpz0(kiv,kjv) = zpzij;
                    et = toc(jTimer);
                    fprintf('Elapsed time is %3.2f seconds\n',et)
                end
                zpz = zpz + zpz0;
            end
            if mv.ecm
                parfor i=1:length(mv.eVecs)
                    fprintf('Starting crossproducts for %20s ... ',mv.eVecs(i).name);
                    iTimer = tic
                    k = 1 + (nTerms + 1)*mv.nPrices + i;
                    zpz0 = zeros(m);
                    zpz0(k,1) = Crossproduct.splxconst(mv.eVecs(i),1,first,last);
                    for jv=1:mv.nPrices
                        j0 = 1 + (jv-1)*(1+nTerms);
                        zpz0(k, j0+1) = Crossproduct.splxspd(mv.eVecs(i),1,mv.dPrices(jv),0,first,last);
                        zpz0(k, j0+1+(1:nTerms)) = Crossproduct.polyxspl(mv.dPrices(jv),mv.polys,1,mv.eVecs(i),1,first,last)';
                    end
                    for j=1:i
                        jj = 1 + (1+nTerms)*mv.nPrices + j;
                        zpz0(k,jj) = Crossproduct.splxspl(mv.eVecs(i),1,mv.eVecs(j),1,first,last);
                    end
                    zpz = zpz + zpz0;
                    et = toc(iTimer);
                    fprintf('Elapsed time is %3.2f seconds\n',et)
                end
            end
            %   fill in upper triangle of zpz.
            zpz = tril(zpz,-1) + diag(diag(zpz)) + tril(zpz,-1)';
            mv.zpz = zpz;
        end
        function [s,t]=zpzLayoutPDL(mv,fDisplay)
            %   The columns of s are
            %   index price# ecm# start end
            if nargin==1; fDisplay=false; end
            nPrices = mv.nPrices;
            nEcm = 0;
            if mv.ecm; nEcm = nPrices-1; end
            s(1).type = 'const';
            s(1).tIndex = 0;
            s(1).first = 1;
            s(1).last = 1;
            for i=1:nPrices
                k = i+1;
                s(k).type = 'price';
                s(k).tIndex = i;
                s(k).first = (i-1)*(1+mv.polys.nTerms) + 2;
                s(k).last = s(k).first + mv.polys.nTerms;
            end
            if mv.ecm
                for iEcm=1:nEcm
                    k = 1 + nPrices + iEcm;
                    s(k).type = 'ecm';
                    s(k).tIndex = iEcm;
                    s(k).first = 1 + nPrices*(1+mv.polys.nTerms) + iEcm;
                    s(k).last = s(k).first;
                end
            end
            if fDisplay
                fprintf('zpzLayoutPDL s:\n')
                disp(struct2table(s))
            end
            k = 0;
            for j=1:length(s);
                for i=j:length(s);
                    k = k + 1;
                    % fprintf('i %d j %d k %d\n',i,j,k)
                    t(k).k=k;
                    t(k).i=i;
                    t(k).j=j;
                end
            end
            if fDisplay
                fprintf('zpzLayoutPDL t:\n')
                disp(struct2table(t))
            end
        end
        function buildzpzSparsePDLpar(mv,mDisplay)
            %   build zpz using calls to sparse matrix routines. 
            %   This calculation uses the zpzLayoutPDL routine to map the cross-product matrix so that most of the matrix (except for the first
            %   column) is computed in one loop. Either 'for' or 'parfor' can be used for this loop.
            if nargin==1; mDisplay=0; end
            if mDisplay>0; fprintf('buildzpzSparsePDLpar\n'); end
            [s,t] = mv.zpzLayoutPDL(false);
            first = mv.firstValid;
            last = mv.lastValid;
            nPrices = mv.nPrices;
            nTerms = mv.polys.nTerms;
            m = length(mv.zNames);
            zpz = zeros(m);
            %   first column (involving constant) is computed directly.
            zfc = NaN(m,1); 
            zfc(1) = last-first+1;
             for i=1:mv.nPrices
                k = (i-1)*(1+mv.polys.nTerms) + 1;
                zfc(k+1) = Crossproduct.spdxconst(mv.dPrices(i),0,first,last);
                kk = k+2 : k+mv.polys.nTerms+1;
                zfc(kk) = Crossproduct.polyxconst(mv.dPrices(i),mv.polys,1,first,last);
            end
            if mv.ecm
                for i=1:length(mv.eVecs)
                    iTimer = tic;
                    k = 1 + (nTerms + 1)*mv.nPrices + i;
                    zfc(k) = Crossproduct.splxconst(mv.eVecs(i),1,first,last);
                end
            end
            zpz(:,1) = zfc;
            dPrices = mv.dPrices.copy();
            polys = mv.polys.copy();
            pTimer = tic;
            pET = 0;
            timingDetailPar = mDisplay;
            parfor k=length(s)+1:length(t)  %   The first column is constants; start with the SECOND column
                jTimer=tic;
                zpz0 = zeros(m);
                i = t(k).i;
                j = t(k).j;
                if mDisplay==2; fprintf('k:%3d i:%d %5s %d j:%d  %5s %d ',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex); end
                switch(s(i).type)
                    case 'price'
                        switch(s(j).type)
                            case 'price'
                                iPrice = s(i).tIndex;
                                jPrice = s(j).tIndex;
                                zpz0(s(i).first,s(j).first) = Crossproduct.spdxspd(dPrices(iPrice), 0, dPrices(jPrice), 0, first, last);
                                zpz0(s(i).first+1:s(i).last, s(j).first) = Crossproduct.polyxspd(dPrices(iPrice), polys, 1, dPrices(jPrice), 0, first,last);
                                zpz0(s(i).first, s(j).first+1:s(j).last) = Crossproduct.polyxspd(dPrices(jPrice), polys, 1, dPrices(iPrice), 0, first,last);
                                zpz0(s(i).first+1:s(i).last, s(j).first+1:s(j).last) = Crossproduct.polyxpoly(dPrices(iPrice), polys, 1, dPrices(jPrice), polys, 1, first, last);
                            case 'ecm'
                                error (sprintf('k:%3d i:%d %5s %d j:%d  %5s %d ',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex));
                            otherwise
                                error (sprintf('k:%3d i:%d %5s %d j:%d  %5s %d ',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex));
                        end
                    case 'ecm'
                        switch(s(j).type)
                            case 'price'
                                iEcm = s(i).tIndex;
                                jPrice = s(j).tIndex;
                                zpz0(s(i).first ,s(j).first) = Crossproduct.splxspd(mv.eVecs(iEcm),1,mv.dPrices(jPrice),0,first,last);
                                zpz0(s(i).first, s(j).first+1:s(j).last) = Crossproduct.polyxspl(mv.dPrices(jPrice),mv.polys,1,mv.eVecs(iEcm),1,first,last)';
                            case 'ecm'
                                iEcm = s(i).tIndex;
                                jEcm = s(j).tIndex;
                                zpz0(s(i).first,s(j).first) = Crossproduct.splxspl(mv.eVecs(iEcm),1,mv.eVecs(jEcm),1,first,last);
                            otherwise
                                error (sprintf('k:%3d i:%d %5s %d j:%d  %5s %d ',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex));
                        end
                    otherwise
                        error (sprintf('k:%3d i:%d %5s %d j:%d  %5s %d ',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex));
                end
                zpz = zpz + zpz0;
                et = toc(jTimer);
                pET = pET + et;
                if mDisplay==2; fprintf('Elapsed time is %3.2f seconds\n',et); end
            end
            if mv.showTiming; fprintf('Within parfor, total elapsed time on all threads is %3.2f sec.\n',pET); end
            pET = toc(pTimer);
            if mv.showTiming; fprintf('Across the full parfor loop, the start to finish elapsed time is %3.2f sec.\n',pET); end
            %   fill in upper triangle of zpz.
            zpz = tril(zpz,-1) + diag(diag(zpz)) + tril(zpz,-1)';
            mv.zpz = zpz;
        end
        function buildzpzSparsePolysPar(mv)
            %   build zpz using calls to sparse matrix routines. Optimized for parallel processing
            first = mv.firstValid;
            last = mv.lastValid;
            m = length(mv.zNames);
            zpz = zeros(m);
            zpz(1,1) = mv.lastValid-mv.firstValid+1;
            %   first column ...
            for i=1:mv.nPrices
                k = (i-1)*(1+mv.polys.nTerms) + 1;
                zpz(k+1,1) = Crossproduct.spdxconst(mv.dPrices(i),0,first,last);
                kk = k+2 : k+mv.polys.nTerms+1;
                zpz(kk,1) = Crossproduct.polyxconst(mv.dPrices(i),mv.polys,1,first,last);
            end
            %   blocks across prices
            nTerms = mv.polys.nTerms;
            dPrices = mv.dPrices.copy;
            polys = mv.polys.copy;
            nPrices = mv.nPrices;
            cdp = parallel.pool.Constant(dPrices);
            cPolys = parallel.pool.Constant(polys);
            np = nPrices*(nPrices+1)/2;
            parfor i=1:np
                iv = floor((1+sqrt(8*i-7))/2)
                jv = i - iv*(iv-1)/2
                k = 1 + (iv-1)*(1+nTerms);
                kiv = (k+1):(k+nTerms+1);
                zpz0 = zeros(m);
                % for jv=1:iv (the loops over rows and columns of the zpz matrix are replaced by a single
                % outer parfor loop.
                tic
                fprintf('Starting iv=%d (%12s) jv=%d (%12s) ... ',iv,cdp.Value(iv).name,jv,cdp.Value(jv).name)
                k = 1 + (jv-1)*(1+nTerms);
                kjv = (k+1):(k+nTerms+1);
                zpzij = zeros(1+nTerms);
                zpzij(1,1) = Crossproduct.spdxspd(cdp.Value(iv), 0, cdp.Value(jv), 0, first, last);
                zpzij(2:end,1) = Crossproduct.polyxspd(cdp.Value(iv), cPolys.Value, 1, cdp.Value(jv), 0, first,last);
                zpzij(1,2:end) = Crossproduct.polyxspd(cdp.Value(jv), cPolys.Value, 1, cdp.Value(iv), 0, first, last);
                zpzij(2:end,2:end) = Crossproduct.polyxpoly(cdp.Value(iv), cPolys.Value, 1, cdp.Value(jv), cPolys.Value, 1, first, last);
                zpz0(kiv,kjv) = zpzij;
                toc
                % end % ... of the jv loop
                zpz = zpz + zpz0;
            end
            if mv.ecm
                eVecs = mv.eVecs;
                cev = parallel.pool.Constant(eVecs);
                necm = length(mv.eVecs);
                parfor i=1:necm
                    zpz0 = zeros(m);
                    fprintf('Starting crossproducts for %s ',cev.Value(i).name);
                    tic
                    k = 1 + (nTerms + 1)*nPrices + i;
                    zpz0(k,1) = Crossproduct.splxconst(cev.Value(i),1,first,last);
                    for jv=1:nPrices
                        j0 = 1 + (jv-1)*(1+nTerms);
                        zpz0(k, j0+1) = Crossproduct.splxspd(cev.Value(i),1,cdp.Value(jv),0,first,last);
                        zpz0(k, j0+1+(1:nTerms)) = Crossproduct.polyxspl(cdp.Value(jv),cPolys.Value,1,cev.Value(i),1,first,last)';
                    end
                    for j=1:i
                        jj = 1 + (1+nTerms)*nPrices + j;
                        zpz0(k,jj) = Crossproduct.splxspl(cev.Value(i),1,cev.Value(j),1,first,last);
                    end
                    zpz = zpz + zpz0;
                    toc
                end
            end
            %   fill in upper triangle of zpz.
            zpz = tril(zpz,-1) + diag(diag(zpz)) + tril(zpz,-1)';
            mv.zpz = zpz;
        end
        function z=buildzpzFull(mv)
            %buildzFull constructs zpz from full (not sparse) data.
            nz = numel(mv.zNames);
            nRows = spd.setgetMax;
            if nRows*nz>1000000; error('buildzpzFull. z matrix has > 1,000,000 entries'); end;
            z = zeros(nRows,nz);
            z(:,1)=1;   %   constant
            k = 1;
            for ip=1:mv.nPrices
                for j=0:mv.arP
                    k = k + 1;
                    z(:,k) = mv.dPrices(ip).lag(j).toCol;
                end
            end
            if mv.ecm
                for i=1:mv.nPrices-1
                    k = k + 1;
                    z(:,k) = mv.eVecs(i).lag.toCol;
                end
            end
            k = ~any(isnan(z),2);
            % fprintf('firstValid=%d lastValid=%d (from full data matrix)\n',find(k,1,'first'),find(k,1,'last'))
            z= z(k,:);
            mv.zpz = z' * z;
            %disptable(z,mv.zNames)
        end
        %% estimation
        function estimateFull(mv)
            % estimateFull: build full data matrices (from sparse) and compute estimates in the usual way.
            z=mv.buildzpzFull;
            %             mv.zpz=z'*z;
            y = z(:,mv.kzlhs);
            %             disp(['Variance of y: ' num2str(var(y))])
            x = z(:,mv.kzrhs);
            xpx = x'*x;
            xpxi = inv(xpx);
            xpy = x'*y;
            b = xpxi*xpy;
            mv.b=b;
            ypy = y'*y;
            e = y-x*b;
            %             disp(['mean of estimated errors:' num2str(mean(e))])
            %             disp('error variance')
            %             disp(cov(e))
            mv.eCov = (ypy - xpy'*b - b'*xpy + b'*xpx*b)/size(y,1);
            v = diag(mv.eCov);
            sim = diag(sqrt(1./v));
            mv.eCorr = sim * mv.eCov * sim;
            mv.seb = sqrt( kron(v', diag(xpxi)) );
            mv.tb = b./mv.seb;
        end
        function estimateSparse(mv)
            % estimate from sparse zpz.
            nObs = mv.zpz(1,1);
            xpx = mv.zpz(mv.kzrhs, mv.kzrhs);
            xpy = mv.zpz(mv.kzrhs, mv.kzlhs);
            xpxi = inv(xpx);
            b = xpxi*xpy;
            mv.b=b;
            ypy = mv.zpz(mv.kzlhs, mv.kzlhs);
            mv.eCov = (ypy - xpy'*b - b'*xpy + b'*xpx*b)/nObs;
            v = diag(mv.eCov);
            sim = diag(sqrt(1./v));
            mv.eCorr = sim * mv.eCov * sim;
            mv.seb = sqrt( kron(v', diag(xpxi)) );
            mv.tb = b./mv.seb;
        end
        %% display, construct and/or format the AR coefficients
        function dispEstimates(mv,eTitle)
            if nargin==1; eTitle=sprintf('VAR/VECM estimates for secondFactor=%d',mv.secondFactor); end
            disp(eTitle)
            t = cellstr(repmat("t  ",1,mv.nPrices));
            colNames = reshape([mv.yNames; t],1,[]);
            %             disptable(mv.b,mv.yNames,mv.xNames);
            disptable(mv.bTable{:,:},colNames,mv.xNames)
            %             disp(mv.bTable);  % use when not sending to diary file
            fprintf('eCov for secondFactor=%d\n',mv.secondFactor);
            disptable(mv.eCov,cellstr(mv.yNames),cellstr(mv.yNames));
            fprintf('eCorr for secondFactor=%d\n',mv.secondFactor);
            disptable(mv.eCorr,cellstr(mv.yNames),cellstr(mv.yNames));
        end
        function b=bTable(mv)
            h={};
            for i=1:mv.nPrices
                h = [h mv.yNames(i) ['t' int2str(i)]];
            end
            btb = reshape([mv.b; mv.tb],[],2*mv.nPrices);
            b=array2table(btb,'VariableNames',cellstr(h),'RowNames',cellstr(mv.xNames));
            % phi = reshape(b(1:2*nLag,:),nLag,2,2);
            % phi = shiftdim(phi,1)
            % gamma = b(end,:)';
        end
        function phi=phiSet(mv,rDisplay) % extract the VAR coefficients
            if nargin<2; rDisplay=false; end
            switch mv.spec
                case Specification.AR
                    phi = zeros(mv.nPrices, mv.nPrices, mv.arP);
                    for ip=1:mv.nPrices
                        for jp=1:mv.nPrices
                            k = (jp-1)*mv.arP + (1:mv.arP);
                            if mv.intercept; k=k+1; end
                            phi(ip,jp,:) = mv.b(k,ip);
                        end
                    end
                case Specification.PDL
                    maxLag = mv.polys.maxLength;
                    phi = zeros(mv.nPrices,mv.nPrices,maxLag);
                    c = mv.polys.designMatrix;
                    nsp = mv.polys.nTerms;
                    for i=1:mv.nPrices      %lhs (equation)
                        for j=1:mv.nPrices  %rhs (variable)
                            ki = (j-1)*nsp + (1:nsp);
                            if mv.intercept; ki = ki + 1; end
                            phi(i,j,:) = c' * mv.b(ki,i);
                        end
                    end
            end
        end      
        function phi=phiSetAR(mv,rDisplay) % extract the VAR coefficients
            %   sets phi to phi(1) phi(2) phi(3) ... where phi(.) is mv.nPrices x mv.nPrices
            if nargin<2; rDisplay=false; end
            phi = zeros(mv.nPrices, mv.nPrices, mv.arP);
            for ip=1:mv.nPrices
                for jp=1:mv.nPrices
                    k = (jp-1)*mv.arP + (1:mv.arP);
                    if mv.intercept; k=k+1; end
                    phi(ip,jp,:) = mv.b(k,ip);
                end
            end
            if rDisplay;
                disp('phi:')
                disp(phi)
            end;
        end
        function phi=phiSetPDL(mv,rDisplay)
            %   Generate the VAR coefficients from polynomial estimates
            if nargin<2; rDisplay=false; end
            maxLag = mv.polys.maxLength;
            phi = zeros(mv.nPrices,mv.nPrices,maxLag);
            c = mv.polys.designMatrix;
            nsp = mv.polys.nTerms;
            %             c = c(mv.kpSelect,:);
            %             nsp = sum(mv.kpSelect); %    The number of poly coefficients actually estimated.
            for i=1:mv.nPrices      %lhs (equation)
                for j=1:mv.nPrices  %rhs (variable)
                    ki = (j-1)*nsp + (1:nsp);
                    if mv.intercept; ki = ki + 1; end
                    phi(i,j,:) = c' * mv.b(ki,i);
                end
            end
            if rDisplay;
                disp('phi:')
                disp(phi)
            end;
            %             mv.phi = phi;
            %             mv.phiLabels = cell(mv.nPrices,mv.nPrices,maxLag);
            %             for i=1:mv.nPrices
            %                 for j=1:mv.nPrices
            %                     for l=1:maxLag
            %                         mv.phiLabels{i,j,l} = sprintf('phi%02d_%02d_%02d',l,i,j);
            %                     end
            %                 end
            %             end
        end
        function vma1=thetaSum(mv) % uses granger representation
            gamma = mv.b(size(mv.b,1)-mv.nPrices+2:end,:)';
            B = [ones(mv.nPrices-1,1) -eye(mv.nPrices-1)];
            gammaNull = null(gamma');
            Bnull = null(B);
            phi = mv.phiSetPDL;
            phiSum = sum(phi,3);
            vma1 = Bnull*inv(gammaNull'*(eye(4)-phiSum)*Bnull)*gammaNull';
        end
        function phiPlot(mv)
            % Plot AR coefficients
            for i=1:mv.nPrices
                figure
                plot(squeeze(mv.phi(i,:,:))')
                title([mv.prices(i).name ' equation']);
                xlabel('Lag')
                ylabel('Coefficient')
                refline(0,0)
                legend(arrayfun(@(x) [x.name ' coeffs'], mv.prices, 'UniformOutput', false))
            end
        end
        %% Simulate
        function simSimple(mv,T) % Simulate a simple Roll model
            spd.setgetMax(T);
            u = randn(1,T);
            p = cumsum(u)+.5*u;
            p1=spl(1:T,p,1,T,'pSim');
            mv.prices = p1;
            mv.nPrices = 1;
        end
        function sim(mv,T,n,nPrices)
            %   T   upper limit on time arrays
            %   n   number of actual price points (density of price matrix is approx n/T)
            %   nPrices
            lambda = n/T;    %   mean arrival rate (proportion of times that will have price changes)
            spd.setgetMax(T);
            rng(8765);
            %   Matlab's exprnd function is parameterized by the mean arrival time; lambda is the arrival rate.
            t=ceil(cumsum(exprnd(1/lambda,1,ceil(2*T*lambda))));    %   simulate some random arrival times
            t=unique(t);    % Select the unique times between 1 and T
            t=t(1<=t & t<=T);
            p=10+cumsum(randn(size(t)));    %   generate prices
            p1 = spl(t,p,t(1),T);
            p1.name='p1';
            mv.prices = spl.empty();
            mv.prices = [mv.prices p1];
            p1.checkForm;
            for j=2:nPrices
                %   Generate some delays as a mixture of exponentials
                t2=t;
                nt = length(t2);
                dLambda=lambda*[1 .1];
                %     dLambda=.2;
                d0=NaN(nt,1);
                for i=nt:-1:1
                    d1 = 2*t2(end);
                    if i<nt; d1 = t2(i+1)-t2(i); end;
                    k = ceil(rand()*length(dLambda));
                    u = rand()*expcdf(d1,1/dLambda(k));
                    d0(i) = expinv(u,1/dLambda(k));
                    % fprintf('i=%d d1=%f d0=%f\n',i,d1,d0);
                    t2(i) = floor(t2(i)+d0(i));
                end
                % [1:nt;t;t2]
                % [mean(d0), min(d0) max(d0)]
                % histogram(d0,'Normalization','probability');
                % hold on;
                % x=.1:0.1:max(d0)
                % p=exppdf(x,mean(d0));
                % plot(x,p,'LineWidth',1.5);
                k = t2<=T; % Trim t2
                px=spl(t2(k),p1.v(k),t2(1),T);
                %px.v = p1.v(k)-j/nPrice;
                %px.v = p1.v(k)-5;
                px.name=['p' int2str(j)];
                mv.prices = [mv.prices px];
            end
            mv.nPrices = nPrices;
        end
        function comparefs(mv)
            %   Verify that full and sparse estimations give the same results.
            mv.buildzpzFull;
            zpz1 = mv.zpz;
            mv.estimateFull
            mv.dispEstimates
            bFull = mv.b;
            nAhead=5;
            irfFull =mv.irfBuildFull(nAhead);
            
            return
            
            mv.buildzpzSparsePolys;
            zpz2 = mv.zpz;
            mv.estimateSparse;
            bSparse = mv.b;
            irfSparse = mv.irfBuildDiff2(nAhead,true);
            
            k = zpz2~=0;
            fprintf(' Largest element of zpz1: %g\n',max(reshape(abs(zpz1(k)),[],1)));
            fprintf('Smallest element of zpz1: %g\n',min(reshape(abs(zpz1(k)),[],1)));
            mad = max(reshape(abs(zpz1(k)-zpz2(k)),[],1));
            [r,c] = find(abs(zpz1-zpz2)==mad & zpz2~=0,1);
            fprintf('\nzpzFull - zpzSparse nmax abs dif: %g at position %d, %d.\n', mad,r,c)
            fprintf('\nbFull-bSparse max abs diff: %g\n',max(abs(reshape(bFull-bSparse,[],1))))
            fprintf('\nirfFull - irfSparse max abs diff: %g\n', max(reshape(abs(irfFull-irfSparse),1,[])))
            
        end
    end
    methods(Static)
        function [irf,irft]=irfCompress(irf0,irft0,irfMax)
            size(irf0)
            size(irft0)
            if nargin<3; irfMax=5000; end
            compressed = any(diff(irft0)~=1);     %   Is the irf already compressed?
            nAhead = irft0(end);
            if compressed | nAhead<=irfMax
                irf = irf0;
                irft = irft0;
            else
                k = unique(reshape((0:100)'*(10.^(0:6)),1,[]));
                k = unique([k nAhead]);
                k = k(k<=nAhead);
                irft = k;
                irf = irf0(:,k+1,:);
            end
        end
    end
    methods(Access=protected)
        function cp=copyElement(o)
            cp = MVARg;
            cp.symbol = o.symbol;
            cp.inDate = o.inDate;
            cp.showTiming = o.showTiming;
            cp.secondFactor = o.secondFactor
            cp.prices = o.prices.copy;
            cp.dPrices = o.dPrices.copy;
            cp.nPrices = o.nPrices;
            cp.spec = o.spec;
            cp.polys = o.polys.copy;
            cp.arP = o.arP;
            cp.intercept = o.intercept;
            cp.ecm = o.ecm;
            cp.firstValid = o.firstValid;
            cp.lastValid = o.lastValid;
            cp.kzlhs = o.kzlhs;
            cp.kzrhs = o.kzrhs;
            cp.priceNames = o.priceNames;
            cp.zNames = o.zNames;
            cp.xNames = o.xNames;
            cp.yNames = o.yNames;
            cp.eVecs = o.eVecs.copy;
            cp.eVecMeans = o.eVecMeans;
            cp.zpz = o.zpz;
            cp.b = o.b;
            cp.bRowNames = o.bRowNames;
            cp.bColNames = o.bColNames;
            cp.seb = o.seb;
            cp.tb = o.tb;
            cp.eCov = o.eCov;
            cp.eCorr = o.eCorr;
            cp.rwd = o.rwd.copy;
            cp.irf = o.irf;
            cp.irft = o.irft;
            cp.irfPacked = o.irfPacked;b
            cp.bIRF = o.bIRF;
            cp.bIRFt = o.bIRFt;
        end
    end
end
