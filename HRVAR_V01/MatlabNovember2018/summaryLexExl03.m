%   Bridged IRFs and related computations for listing/non-listing exchanges
%
clc; 
clf; 
clear all
close all force
format; format compact
addpath('./mClasses','./mFiles')

fontSize=9;
shockVar='nbbLex';
impactVar='nbbExl';
inDir = '.\LexExl03\save';
fileNameFmt = 'MVARiLexExl03_%s_20161003.mat';
diary off;
diaryFile = 'summaryLexExl03.txt';
delete(diaryFile)
diary(diaryFile)
symbols ={'IBM','NVDA'};
for iSymbol=1:2
    symbol = symbols{iSymbol};
    inFile = sprintf(fileNameFmt,symbol);
    inFileFull = fullfile(inDir,inFile);
    fprintf('\n%s\nResults for %s\n',repmat('=',1,50),inFileFull)
    load(inFileFull)
    %%   Plot individual IRFs
    %%   Construct bridged IRF
    % mva.irfPlot('nbbPart','nbbSIP')
    kBridge = 10:10:100;   %   bridging indexes
    %   bridged impulse response functions (over one time scale)
    mva = flip(mva);    %   reverse so that mva is ordered from highest resolution to lowest
    mv = mva(1);
    for i=1:length(mva)
        fprintf('\n%s\n',repmat('-',1,30))
        mva(i).rwd.rDisplay
    end
    nShocks = mv.nPrices;
    nShockDisplay = 1;
    k = mv.irft<=kBridge(end);
    mv.bIRF = mv.irfPacked(:,k,:);
    mv.bIRFt = mv.irft(k);
    colNames = strsplit(int2str(mv.bIRFt));
    fprintf('%s Bridged IRF starting at secondFactor=%-6d\n',symbol,mv.secondFactor)
    % disptable(squeeze(mv.bIRF(:,:,1)),colNames,mv.priceNames)
    %   move to coarser time scales
    for iSec=2:length(mva)
        mv1 = mva(iSec-1);
        mv = mva(iSec);
        if mv.secondFactor==1; break; end
        fprintf('%s\nISEC=%d secondFactor=%d\n',repmat('-',1,50),iSec,mv.secondFactor)
        %         colNamesLast = strsplit(int2str(b(iSec-1).irft));
        nAhead = kBridge(end);
        if mv.secondFactor==10; nAhead = 500*mv.secondFactor; end
        phi = mv.phiSet;
        irf = zeros(mv.nPrices,nAhead+1,mv.nPrices);
        for iShock=1:mv.nPrices
            % fprintf('bridged IRF ... iShock=%d\n',iShock)
            [la, lb]=ismember(kBridge, mv1.bIRFt);
            if any(la==0); error('iShock=%d. Can''t find elements in mv1.bIRFt',iShock); end
            pLag = mv1.bIRF(:,[1 lb],iShock);
            f = mv.vForecast(pLag,phi,nAhead);
            f = [pLag pLag(:,end)+cumsum(f,2)];
            f = f(:,1:nAhead+1);
            irf(:,:,iShock) = f;
        end
        irft = 0:nAhead;
        [la,lb] = ismember(mv.irft,irft);
        mv.bIRFt = irft(lb(la));
        mv.bIRF = irf(:,lb(la),:);
        %   Display sequence of irfs
        % for iShock=1:nShockDisplay
        %     fprintf('Original irf:\n')
        %     colNames=strsplit(int2str(mv.irft));
        %     m = squeeze( mv.irfPacked(:,:,iShock) );
        %     disptable(m,colNames,mv.priceNames)
        %     fprintf('Bridged irf:\n')
        %     s2=strsplit(int2str(mv.bIRFt));
        %     m2 = squeeze( mv.bIRF(:,:,iShock) );
        %     disptable(m2,s2,mv.priceNames)
        % end
    end
    fprintf('\n%s\nBridged random walk decomposition for %s\n',repmat('-',1,50),symbol)
    r = randomWalkDecomp;
    bSecondFactors = arrayfun(@(x) x.secondFactor, mva);
    kCoarse = find(10==bSecondFactors);
    if isempty(kCoarse); error('secondFactor=10 not found in %s',int2str(bSecondFactors)); end
    r.init(mva(1));
    mvCoarse = mva(kCoarse);
    vma1 = mean(squeeze(mvCoarse.bIRF(:,end,:)),1)
    r.solve(vma1)
    r.isBoundsGrouped(mva(1).rwd.isGroups)
    r.rDisplay
    %% Plot bridged IRFs
    mv = mva(1);
    [~,iShock]=ismember(shockVar,mv.priceNames);
    [~,jVariable]=ismember(impactVar,mv.priceNames);
    irfTitle = sprintf('%s bridged IRF: %s > %s\n',symbol,shockVar,impactVar);
    fontSize=10;
    figure('Units','inches','Position',[1,1,6,4]);
    plotSet1(irfTitle,fontSize)
    hold
    yMax = 0;
    for iSec=1:kCoarse
        mv = mva(iSec);
        d = (iSec-1)*.01;
        x = mv.bIRFt(2:end)/mv.secondFactor;
        y = d+squeeze(mv.bIRF(jVariable,2:end,iShock));
        semilogx(x,y);
        yMax = max([yMax y]);
        ax = gca;
        % hold on;
        % switch mv.secondFactor
        %     case 1; pText='1 sec';
        %     case 10; pText='100 ms';
        %     case 100; pText='10 ms';
        %     case 1000; pText='1 ms';
        %     case 10000; pText='100 \mus';
        %     case 100000; pText='10 \mus';
        %     otherwise; pText=''
        % end
        % dx = 0;
        % dy = -(ax.YLim(2)-ax.YLim(1))/30;
        % X = mv.bIRFt(2)/mv.secondFactor + dx;
        % Y = mv.bIRF(jVariable,2,iShock) + dy;
        % text(X,Y,[pText],'FontSize',fontSize);
    end
    
    ytix = arrayfun(@(x) sprintf('$%4.2f',x),ax.YTick,'UniformOutput',false);
    k=ax.YTick<0;
    if any(k); ytix{k}=''; end
    ax.YTickLabel = ytix;
    ax.YLim(1) = min(ax.YLim(1),-0.01);
    ax.YLim(2) = yMax+.05;
    ax.YTickMode='manual';  %   'manual' prevents automatic redrawing of ticks if the graph is resized.
    lineSet1
    hold off;   
    gFile = sprintf('summaryLexExl03 %4s',symbol);
    saveas(gcf,gFile,'fig')
    print('-r300',gFile,'-djpeg')
end
diary off