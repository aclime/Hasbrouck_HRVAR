classdef randomWalkDecomp < matlab.mixin.Copyable
    %randomWalkDecomp constructs a random-walk decomposition
    properties
        nPrices=0
        secondFactor=[] %   used in annualization of sdw
        vma1=[]
        eCov=[]
        yNames={}
        var_w=[]
        var_wAnnual=[]
        isBounds=[]
        allPerms=[]
        vwAllPerms=[]
        isGroups=[]
        isGroupBounds=[]
        isGroupLabels=[]
    end
    methods
        function init(o, mv)    %   initialize from an MVARi object
            if ~isa(mv,'MVARi'); error('init expecting to be called with MVARi object'); end
            o.nPrices = mv.nPrices;
            o.secondFactor = mv.secondFactor;
            o.yNames = mv.priceNames;
            o.eCov = mv.eCov;
            vma1 = mv.thetaSum;
            o.vma1 = vma1(1,:);
            o.var_w = o.vma1 * o.eCov * o.vma1'; %   For a price discovery analysis, all entries should be equal.
            o.var_wAnnual = o.var_w*(250*6.5*3600*mv.secondFactor);
            o.buildAllPerms;
            is = o.vwAllPerms/mean(sum(o.vwAllPerms,2));
            o.isBounds = vertcat(min(is,[],1),max(is,[],1));
        end
        %   Solve with alternate values of vma1 and eCov
        %   (typically used when vma1 is constructed from a multiscale irf)
        function solve(o,vma1,eCov)
            if nargin>=2 && ~isempty(vma1); o.vma1 = vma1; end
            if nargin>=3 && ~isempty(eCov); o.eCov = eCov; end
            o.var_w = o.vma1 * o.eCov * o.vma1'; %   For a price discovery analysis, all entries should be equal.
            if ~isempty(o.secondFactor) o.var_wAnnual = o.var_w*(250*6.5*3600*o.secondFactor); end
            o.buildAllPerms;
            is = o.vwAllPerms/mean(sum(o.vwAllPerms,2));
            o.isBounds = vertcat(min(is,[],1),max(is,[],1));
        end
        function initm(o,vecm,secondFactor)  %   initialize from a matlab vecm object
            if ~isa(vecm,'vecm'); error('Expected an (estimated) toolbox vecm object.'); end
            if isempty(vecm.Adjustment); error('toolbox vecm (passed argument) has not been estimated.'); end
            if nargin>=3; o.secondFactor=secondFactor; end
            o.nPrices = vecm.NumSeries;
            o.yNames = cellstr(vecm.SeriesNames);
            o.eCov = vecm.Covariance;
            %   From granger representation:
            gamma = vecm.Adjustment;
            gammaNull = null(gamma');
            B = vecm.Cointegration';
            Bnull = null(B);
            phi = reshape(cell2mat(vecm.ShortRun),o.nPrices,o.nPrices,[]);
            phiSum = sum(phi,3);
            vma1 = Bnull*inv(gammaNull'*(eye(o.nPrices)-phiSum)*Bnull)*gammaNull';
            o.vma1 = vma1(1,:);
            o.var_w = o.vma1 * o.eCov * o.vma1'; %   For a price discovery analysis, all entries should be equal.
            if ~isempty(o.secondFactor); o.var_wAnnual = o.var_w*(250*6.5*3600*secondFactor); end
            o.buildAllPerms;
            is = o.vwAllPerms/mean(sum(o.vwAllPerms,2));
            o.isBounds = vertcat(min(is,[],1),max(is,[],1));
        end
        function rDisplay(r,title)
            if nargin==1    % no title; make one up.
                if isempty(r.secondFactor);
                    title="Random walk decomposition in event time";
                else
                    title=sprintf('Random walk decomposition for secondFactor=%d',r.secondFactor);
                end
            end
            disp(title);
            disp('Sum of vma coefficients:')
            disptable(r.vma1(1,:),r.yNames);
            fprintf('per period var_w:%12.6f sd_w:%12.6f\n',r.var_w,sqrt(r.var_w));
            if ~isempty(r.var_wAnnual)
                fprintf('annualized var_w:%12.6f sd_w:%12.6f\n',r.var_wAnnual,sqrt(r.var_wAnnual));
            end
            if ~isempty(r.isBounds)
                disp('info share bounds:');
                disptable(r.isBounds',{'Min','Max'},r.yNames,'%8.5f');
            end
            if ~isempty(r.isGroupBounds)
                disp('Grouped info share bounds')
                disptable(r.isGroupBounds,{'Min','Max'},r.isGroupLabels,'%8.5f'); end
        end
        function isBoundsGrouped(r,groups)
            is = r.vwAllPerms/mean(sum(r.vwAllPerms,2));
            nGroups = length(groups);
            groupLabels = cell(1,2);
            isb = NaN(nGroups,2);
            for i=1:nGroups
                [k,ll] = ismember(groups{i},r.yNames);
                if any(~k)
                    disp(groups{i})
                    disp('not found in')
                    disp(r.yNames)
                    error('isBoundsGrouped')
                end
                isg = sum(is(:,ll(k)),2);
                isb(i,:) = [min(isg,[],1) max(isg,[],1)];
                groupLabels{i} = char(join(groups{i}));
            end
            r.isGroups = groups;
            r.isGroupBounds = isb;
            r.isGroupLabels = groupLabels;
        end
        function [vw, infoShare]=rwPerm(r,jPerm)
            vma = mean(r.vma1,1);
            eCovPerm = r.eCov(jPerm,jPerm);
            vmaPerm = vma(jPerm);
            L = chol(eCovPerm,'lower')
            llp = L*L'
            lpl = L'*L
            w = (vmaPerm*L);
            vw = w.^2;
            infoShare = vw/sum(vw);
        end
        function buildAllPerms(r)
            vma = mean(r.vma1,1);
            p = perms(1:r.nPrices);
            r.vwAllPerms = NaN(size(p));
            nPerms = size(p,1);
            for ip=1:nPerms
                k = p(ip,:);
                eCovPerm = r.eCov(k,k);
                vmaPerm = vma(k);
                L = chol(eCovPerm,'lower');
                w = vmaPerm*L;
                r.vwAllPerms(ip,k) = w.^2;
            end
            r.allPerms = p;
        end
    end
    methods(Access=protected)
        function r=copyElement(o)
        r = randomWalkDecomp;
        r.nPrices = o.nPrices;
        r.secondFactor = o.secondFactor;
        r.vma1 = o.vma1;
        r.eCov = o.eCov;
        r.yNames = o.yNames;
        r.var_w = o.var_w;
        r.var_wAnnual = o.var_wAnnual;
        r.infoShareBounds = o.infoShareBounds;
        r.allPerms = o.allPerms;
        r.vwAllPerms = o.vwAllPerms;
        r.isGroups = o.isGroups;
        r.isGroupBounds = o.isGroupBounds;
        r.isGroupLabels = o.isGroupLabels;
        end
    end
end

