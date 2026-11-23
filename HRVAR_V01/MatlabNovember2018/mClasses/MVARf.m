classdef MVARf < matlab.mixin.Copyable
    % Microstructure VAR/VECM class
    % Coefficient structure is standard (non-polynomial)
    properties
        secondFactor = 1        % The number of time units per second.
        prices=spl.empty
        dPrices=spd.empty
        nPrices
        %   Specification variables
        intercept=true  %   include intercept in specifications
        ecm=true        %   include error correction terms
        firstValid %   first valid time (index) in the 'data' matrix
        lastValid  %   last         
        arP;                % order of VAR
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
        seb    %   se's of coefficients
        tb     %   t's of coefficients
        eCov   %   residual covariance matrix
        eCorr  %   residual correlation matrix
        irft    %   time indices for irf array
        irf     %   nPrices x nAhead x nPrices(=nShocks)
        rwd     %   random walk decomposition (object)
    end
    methods
        function o=MVARf(n)
            if nargin~=0; o(n)=MVARf; end
        end
                
        function setup(mv)
            %   take price differences, form cointegrating vectors, etc.
            mv.nPrices = length(mv.prices);
            mv.dPrices(mv.nPrices) = spd;
            for i=1:mv.nPrices
                if mv.prices(i).maxSize~=mv.prices(1).maxSize; error('MVARf.setup. price vectors don''t have the same maxSize.'); end
                mv.dPrices(i) = mv.prices(i).diff;
            end
            if mv.ecm && mv.nPrices>=2
%                 fprintf('Computing eVecs\n')
                mv.eVecs(mv.nPrices-1)=spl;
                for i=1:mv.nPrices-1
%                     fprintf('%s [minus]\n',mv.prices(1).toString)
%                     fprintf('%s [equals]\n',mv.prices(i+1).toString)
                    mv.eVecs(i) = mv.prices(1).splMinus(mv.prices(i+1));
                    mv.eVecs(i).name = [mv.prices(1).name '-' mv.prices(i+1).name];
%                     fprintf('%s\n',mv.eVecs(i).toString)
                end
            end
        end
        function setNames(mv,nDisplay)   %   set up row and column labels for various matrices;
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
        function setValid(mv,maxLag,display)
            if nargin<3; display=false; end
            first=0;
            maxSize = mv.prices(1).maxSize;
            spd.setgetMax(maxSize);
            last = maxSize+1;
            %last=spd.setgetMax+1;
            maxLag = max(1,mv.arP);
            for i=1:mv.nPrices
                if display; fprintf('MVARf.setValid %12s first=%d:',mv.prices(i).name, mv.prices(i).firstValid); end
                if mv.prices(i).maxSize~=maxSize; error('MVARf.setValid. price vectors don''t have the same maxSize.'); end
                first=max([first,mv.prices(i).firstValid+maxLag+1]);
                last=min([last,mv.prices(i).lastValid]);
                if display; fprintf('first=%d last=%d\n',first,last); end
            end
            mv.firstValid = first;
            mv.lastValid = last;
        end
        function eVecDemean(mv)
            if ~mv.ecm; return; end
            mv.eVecMeans = NaN(1,mv.nPrices-1);
            n = mv.lastValid - mv.firstValid + 1;
            for i=1:mv.nPrices-1
                z = mv.eVecs(i).sum(mv.firstValid, mv.lastValid, 1);
                mv.eVecMeans(i) = z/n;
                mv.eVecs(i).v = mv.eVecs(i).v - mv.eVecMeans(i);
            end
        end
        function irfBuild2(mv,nAhead,display) %  Direct IRF computation
            if nargin<3; display=false; end
            mv.phiSet
            maxLag = size(mv.phi,3);    % This may be <= the maxLag in the zpz matrix.
            irf = zeros(mv.nPrices,nAhead,mv.nPrices);
            T = maxLag + nAhead;
            dpf = zeros(mv.nPrices,nAhead);
            mPhi = reshape(mv.phi,mv.nPrices,[]);   %   Restack phi into an nPrice x (nPrice*maxLag) maxtrix
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
                dpf(:,:,iShock) = flip(dp(:,1:nAhead),2);
                irf(:,:,iShock) = cumsum( flip(dp(:,1:nAhead),2), 2);
            end
            mv.irf = irf;
            if display
                nPrint = min(10,nAhead);
                disp('dpIRF:')
                for i=1:mv.nPrices; disptable(dpf(:,1:nPrint,i),strsplit(int2str(0:(nPrint-1))),mv.priceNames); end
                disp('pIRF:')
                for i=1:mv.nPrices; disptable(irf(:,1:nPrint,i),strsplit(int2str(0:(nPrint-1))),mv.priceNames); end
            end
        end
        function irfBuild(mv,nAhead,display) %  Direct IRF computation
            if nargin<3; display=false; end
            if mv.arP>0
                phi=flip(mv.phiSet,3);
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
                if display
                    nPrint = min(nAhead,10);
                    fprintf('MVARf.irfBuild iShock=%d\n',iShock)
                    disptable(dp(:,1:10),strsplit(int2str(1:10)),mv.priceNames);
                end;
                irf(:,:,iShock) = cumsum(dp(:,1+mv.arP:end),2);
            end
            mv.irf = irf;
            mv.irft = 0:size(irf,2)-1;
            if display
                for i=1:mv.nPrices; 
                    x = [irf(:,1:nPrint,i) irf(:,end,i)];
                    s = strsplit(int2str([0:nPrint-1 size(irf,2)-1]));
                    disptable(x,s,mv.priceNames); 
                end
            end
        end
        function buildzpzSparse(mv) %   January, 2018
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
        function buildzpzSparsePolys(mv)
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
            for iv=1:nPrices
                k = 1 + (iv-1)*(1+nTerms);
                kiv = (k+1):(k+nTerms+1);
                zpz0 = zeros(m);
                for jv=1:iv
                    tic
                    fprintf('Starting iv=%d (%12s) jv=%d (%12s) ... ',iv,dPrices(iv).name,jv,dPrices(jv).name)
                    k = 1 + (jv-1)*(1+nTerms);
                    kjv = (k+1):(k+nTerms+1);
                    zpzij = zeros(1+nTerms);
                    zpzij(1,1) = Crossproduct.spdxspd(dPrices(iv), 0, dPrices(jv), 0, first, last);
                    zpzij(2:end,1) = Crossproduct.polyxspd(dPrices(iv), polys, 1, dPrices(jv), 0, first,last);
                    zpzij(1,2:end) = Crossproduct.polyxspd(dPrices(jv), polys, 1, dPrices(iv), 0, first, last);
                    zpzij(2:end,2:end) = Crossproduct.polyxpoly(dPrices(iv), polys, 1, dPrices(jv), polys, 1, first, last);
                    zpz0(kiv,kjv) = zpzij;
                    toc
                end
                zpz = zpz + zpz0;
            end
            if mv.ecm
                for i=1:length(mv.eVecs)
                    fprintf('Starting crossproducts for %s ',mv.eVecs(i).name);
                    tic
                    k = 1 + (nTerms + 1)*mv.nPrices + i;
                    zpz(k,1) = Crossproduct.splxconst(mv.eVecs(i),1,first,last);
                    for jv=1:mv.nPrices
                        j0 = 1 + (jv-1)*(1+nTerms);
                        zpz(k, j0+1) = Crossproduct.splxspd(mv.eVecs(i),1,mv.dPrices(jv),0,first,last);
                        zpz(k, j0+1+(1:nTerms)) = Crossproduct.polyxspl(mv.dPrices(jv),mv.polys,1,mv.eVecs(i),1,first,last)';
                    end
                    for j=1:i
                        jj = 1 + (1+nTerms)*mv.nPrices + j;
                        zpz(k,jj) = Crossproduct.splxspl(mv.eVecs(i),1,mv.eVecs(j),1,first,last);
                    end
                    toc
                end
            end
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
        function dispEstimates(mv,eTitle)
            if nargin==1; eTitle=sprintf('VAR/VECM estimates for secondFactor=%d arP=%d',mv.secondFactor,mv.arP); end
            disp(eTitle)
            % disptable(mv.b,mv.yNames,mv.xNames);
            disp(mv.bTable);
            fprintf('eCov for secondFactor=%d arP=%d\n',mv.secondFactor,mv.arP);
            disptable(mv.eCov,cellstr(mv.yNames),cellstr(mv.yNames));
            fprintf('eCorr for secondFactor=%d arP=%d\n',mv.secondFactor,mv.arP);
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
        function phi=phiSet(mv,display) % extract the VAR coefficients
            %   sets mv.phi to phi(1) phi(2) phi(3) ... where phi(.) is mv.nPrices x mv.nPrices
            if nargin<2; display=false; end
            phi = zeros(mv.nPrices, mv.nPrices, mv.arP);
            for ip=1:mv.nPrices
                for jp=1:mv.nPrices
                    k = (jp-1)*mv.arP + (1:mv.arP);
                    if mv.intercept; k=k+1; end
                    phi(ip,jp,:) = mv.b(k,ip);
                end
            end
            if display;
                disp('phi:')
                disp(mv.phi)
            end;
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
            cp = MVARe;
            cp.secondFactor = o.secondfFactor
            cp.prices = o.prices.copy;
            cp.dPrices = o.dPrices.copy;
            cp.nPrices = o.nPrices;
            cp.polys = o.polys.copy;
            cp.kpSelect = o.kpSelect;
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
%             cp.phi = o.phi;
%             cp.phiLabels = o.phiLabels;
            cp.rwd = o.rwd.copy;
            cp.irf = o.irf;
            cp.irft = o.irft;
        end
    end
end
