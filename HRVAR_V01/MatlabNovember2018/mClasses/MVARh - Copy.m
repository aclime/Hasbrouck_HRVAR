classdef MVARh < matlab.mixin.Copyable
    % Microstructure VAR/VECM class
    % The notional data matrix is z, and the notional crossproduct matrix is z'p (zpz).
    % The column names of z are set in setnamesPDL
    % [1 (P1 block) (P2 block) ... (PnPrice block) eVecs]
    % where a price block includes the price (not lagged) and all polynomials or AR coefficients.
    % the eVecs are the error correction terms.
    properties
        symbol=char.empty  %   ticker symbol
        inDate=[]   %   sample date
        showTiming = true
        secondFactor = 1        % The number of time units per second.
        prices=spl.empty
        dPrices=spd.empty
        nPrices=[]
        %   Specification variables
        spec=Specification.PDL;         %   specification type: AR or PDL. NOTE: AR is not fully implemented.
        intercept=true  %   include intercept in specifications
        ecm=true        %   include error correction terms
        firstValid=[]      %   first valid time (index) in the 'data' matrix
        lastValid=[]       %   last ...
        %   if spec=AR, then arP must be set; if spec=PDL then polys must be set.
        arP=[];                    % order of VAR (NOT yet implemented)
        polys=polynom.empty     % array of polynomials that define the design matrix for the coefficients; set by user.
        kzlhs=[]       %   logical index vector within z to left hand side vars
        kzrhs=[]       %   ... and right hand side vars
        priceNames={}  %   names of price variables (levels, like 'bid' and 'ask'). Set by user.
        zNames={}      %   Variable names corresponding in the cross-product matrix; built in setNames routines
        xNames={}      %   rhs variables in the cross-product matrix; built in setNames routines.
        yNames={}      %   lhs variables in the cross-product matrix; built in setNames routines.
        eVecs=spl.empty   %   error (cointegration) vectors
        eVecMeans=[]   % means of error correction terms
        zpz=[]     %   full cross-product matrix
        %   parameters
        b=[]       %   coefficients
        bRowNames={}
        bColNames={}
        seb=[]     %   se's of coefficients
        tb=[]      %   t's of coefficients
        eCov=[]    %   residual covariance matrix
        eCorr=[]   %   residual correlation matrix
        irft=[]    %   time indices for irf array (intervals; to get to real time, divide by secondFactor)
        irf=[]     %   nPrices x nAhead x nPrices(=nShocks)
        irfPacked=[]   %   Condensed
        bIRFt=[]   %   time indices for bridged irf (intervals; to get to real time, divide by secondFactor)
        bIRF=[]    %   bridged IRF 
        rwd=randomWalkDecomp.empty      %   random walk decomposition in clock time
        %%  Event time estimates
        vecmEst=vecm.empty              %   estimated vecm object from econometric toolbox.
        vecmSE=[]
        nObsEventTime=[]
        rwdet=randomWalkDecomp.empty    %   random walk decomposition in event time.
    end
    methods
        %% initialization and setup routines.
        function o=MVARh(n) % array constructor
            %   MVARh(n) returns an array of MVARh's.
            if nargin~=0; o(n)=MVARh; end
        end
        function setup(mv,mDisplay)
            %   take price differences, form cointegrating vectors, etc.
            if nargin<2; mDisplay=false; end
            if mDisplay; fprintf('setup.\n'); end
            mv.nPrices = length(mv.prices);
            mv.dPrices(mv.nPrices) = spd;
            for i=1:mv.nPrices
                if mv.prices(i).maxSize~=mv.prices(1).maxSize; error('MVARh.setup. price vectors don''t have the same maxSize.'); end
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
            for i=1:mv.nPrices
                n = length(mv.prices(i).i);
                nValid = sum(mv.prices(i).i>=mv.firstValid & mv.prices(i).i<=mv.lastValid);
                fprintf('%12s %10d total observations, %10d in valid range.\n',mv.prices(i).name,n,nValid);
            end
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
        %   'Direct': The irf is computed by iterating the full AR representation, derived from the PDLs. (The alternative is via
        %       the state representation of the cumulative PDLs, described in the computational appendix.)
        %   'Packed': The returned irf is based on a log time scale
        %   Some of the routines use parfor loops (parallel processing, for computational efficiency)
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
        function [irfPacked,irft]=irfPDLpacked(mv,nAhead,displayLevel,nCPU)     % (UPDATED 7/18/18) 
            if nargin<4; nCPU=0; end
            if nargin<3; displayLevel=0; end
            if displayLevel>0; fprintf('irfPDLpacked (nCPU=%d)\n',nCPU); end
            m = ceil(log10(nAhead));
            k = unique( horzcat(reshape( (0:10)'*(10.^(0:m)), 1, []), nAhead) );
            irft = k(k<=nAhead);
%             if fDisplay
%                 fprintf('irfPDLpacked. irft:\n')
%                 disptable(irft,split(int2str(1:length(irft))))
%             end
            irfPacked = NaN(mv.nPrices,length(irft),mv.nPrices);
            fDisplayPar = displayLevel>=4;
            pTimer = tic;
            pET = 0;
            parfor (iShock=1:mv.nPrices, nCPU)
                jTimer = tic;
                e0 = zeros(mv.nPrices,1);
                e0(iShock)=1;
                f = fPDLpacked(mv,e0,nAhead,irft,fDisplayPar);
                irfPacked(:,:,iShock) = f;
                et = toc(jTimer);
                pET = pET+et;
                if displayLevel>=2; fprintf('IRF computation for shock %d. Elapsed time is %3.2f seconds.\n',iShock,et); end
            end
            if displayLevel>=2;
                fprintf('Within parfor, total elapsed time on all threads is %3.2f seconds.\n',pET);
                pET = toc(pTimer);
                fprintf('Across the full parfor loop, the start to finish elapsed time is %3.2f sec.\n',pET);
            end
            if displayLevel>=3
                disp('irfPDLpacked irfPacked');
                for i=1:mv.nPrices
                    disptable(irfPacked(:,:,i),split(int2str(irft)),mv.priceNames);
                end
            end
            mv.irfPacked = irfPacked;
            mv.irft = irft;
            mv.irf = [];
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
        function buildzpzSparsePDL0(mv, nCPU)
            %   build zpz using direct calls to sparse matrix routines (does not use s and t; calls to parfor are not efficiently organized:
            %   parfor loops over rows, and the inner loops are over columns.
            %   nCPU (optional) is used in parallel processing. Set to zero (the default) to use serial execution
            if nargin<2; nCPU=0; end
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
            parfor (iv=1:nPrices, nCPU)
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
                parfor (i=1:length(mv.eVecs),nCPU)
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
        function [s,t]=zpzLayoutPDL(mv,fDisplay)    %   data structures to map linear index into the 2D cross-product matrix.
            %   parfor (parallel looping) does not nest: we can't have a loop over rows nested inside of a loop over columns or vice versa.
            %   t maps a linear index into rows/columns of the crossproduct matrix zpz
            %   s describes the arrangement of variables along a row (or column)
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
                    t(k).i=i;   % In building the zpz matrix, i and j will index into s.
                    t(k).j=j;
                end
            end
            if fDisplay
                fprintf('zpzLayoutPDL t:\n')
                disp(struct2table(t))
            end
        end
        function buildzpzSparsePDLpar(mv,displayLevel,nCPU)
            %   build zpz using calls to sparse matrix routines. 
            %   This calculation uses the zpzLayoutPDL routine to map the cross-product matrix so that most of the matrix (except for the first
            %   column) is computed in one loop.
            %   nCPU controls the parallel processing. Set to zero to force serial execution.
            %   displayLevel: 0 (message printed at invocation only); 2 (print timing details)
            if nargin<3; nCPU=0; end;   %   defaults to serial execution
            if nargin<2; displayLevel=0; end
            if displayLevel>0; fprintf('buildzpzSparsePDLpar (nCPU=%d)\n',nCPU); end
            layoutDisplay = displayLevel>=2;
            [s,t] = mv.zpzLayoutPDL(layoutDisplay);
            first = mv.firstValid;
            last = mv.lastValid;
            nPrices = mv.nPrices;
            nTerms = mv.polys.nTerms;
            m = length(mv.zNames);
            zpz = zeros(m);
            %   first column (crossproduct of each term with the constant) is computed directly.
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
            pTimer = tic;
            pET = 0;
            %ticBytes(gcp)
            parfor (k=length(s)+1:length(t),nCPU)  %   In this loop, the first column is constants; start with the SECOND column
%                   for k=length(s)+1:length(t) %   Use for straight-through looping
%                   for k=27:27 %   Use for checking a particular block
                jTimer=tic;
                zpz0 = zeros(m);
                i = t(k).i;
                j = t(k).j;
                if displayLevel>=2; fprintf('Starting k:%3d i:%d %5s %d j:%d  %5s %d \n',k,i,s(i).type,s(i).tIndex,j,s(j).type,s(j).tIndex); end
                switch(s(i).type)
                    case 'price'
                        switch(s(j).type)
                            case 'price'
                                iPrice = s(i).tIndex;
                                jPrice = s(j).tIndex;
                                if displayLevel>=4; fprintf('k:%3d calling spdxspd.\n',k); end
                                zpz0(s(i).first,s(j).first) = Crossproduct.spdxspd(mv.dPrices(iPrice), 0, mv.dPrices(jPrice), 0, first, last);
                                if displayLevel>=4; fprintf('k:%3d calling polyxspd(i,j)\n',k); end
                                zpz0(s(i).first+1:s(i).last, s(j).first) = Crossproduct.polyxspd(mv.dPrices(iPrice), mv.polys, 1, mv.dPrices(jPrice), 0, first,last);
                                if displayLevel>=4; fprintf('k:%3d calling polyxspd(j,i)\n',k); end
                                zpz0(s(i).first, s(j).first+1:s(j).last) = Crossproduct.polyxspd(mv.dPrices(jPrice), mv.polys, 1, mv.dPrices(iPrice), 0, first,last);
                                if displayLevel>=4; fprintf('k:%3d calling polyxpoly\n',k); end
                                zpz0(s(i).first+1:s(i).last, s(j).first+1:s(j).last) = Crossproduct.polyxpoly(mv.dPrices(iPrice), mv.polys, 1, mv.dPrices(jPrice), mv.polys, 1, first, last);
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
                                if displayLevel>=4; fprintf('k:%3d calling splxspd\n',k); end
                                zpz0(s(i).first ,s(j).first) = Crossproduct.splxspd(mv.eVecs(iEcm),1,mv.dPrices(jPrice),0,first,last);
                                if displayLevel>=4; fprintf('k:%3d calling polyxspl\n',k); end
                                zpz0(s(i).first, s(j).first+1:s(j).last) = Crossproduct.polyxspl(mv.dPrices(jPrice),mv.polys,1,mv.eVecs(iEcm),1,first,last)';
                            case 'ecm'
                                iEcm = s(i).tIndex;
                                jEcm = s(j).tIndex;
                                if displayLevel>=4; fprintf('k:%3d calling splxspl\n',k); end
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
                if displayLevel>=2; fprintf('Ending k:%3d Elapsed time is %3.2f seconds\n',k,et); end
            end
            %tocBytes(gcp)
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
            vma1 = Bnull*inv(gammaNull'*(eye(mv.nPrices)-phiSum)*Bnull)*gammaNull';
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
    methods(Access=protected)
        function cp=copyElement(o)
            cp = MVARh;
            
            cp.prices = o.prices.copy;
            cp.dPrices = o.dPrices.copy;
            cp.polys = o.polys.copy;
            cp.eVecs = o.eVecs.copy;
            cp.rwd = o.rwd.copy;
            
            cp.symbol = o.symbol;
            cp.inDate = o.inDate;
            cp.showTiming = o.showTiming;
            cp.secondFactor = o.secondFactor;
            cp.nPrices = o.nPrices;
            cp.spec = o.spec;
            cp.intercept = o.intercept;
            cp.ecm = o.ecm;
            cp.firstValid = o.firstValid;
            cp.lastValid = o.lastValid;
            cp.arP = o.arP;
            cp.kzlhs = o.kzlhs;
            cp.kzrhs = o.kzrhs;
            cp.priceNames = o.priceNames;
            cp.zNames = o.zNames;
            cp.xNames = o.xNames;
            cp.yNames = o.yNames;
            cp.eVecMeans = o.eVecMeans;
            cp.zpz = o.zpz;
            cp.b = o.b;
            cp.bRowNames = o.bRowNames;
            cp.bColNames = o.bColNames;
            cp.seb = o.seb;
            cp.tb = o.tb;
            cp.eCov = o.eCov;
            cp.eCorr = o.eCorr;
            cp.irft = o.irft;
            cp.irf = o.irf;
            cp.irfPacked = o.irfPacked;
            cp.bIRFt = o.bIRFt;
            cp.bIRF = o.bIRF;
            
        end
    end
end
