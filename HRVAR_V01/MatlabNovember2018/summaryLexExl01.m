%   Listing and non-listing exchanges; summary and graphing
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
inDir = '.\HPC results\LexExl01\save';
fileNameFmt = 'MVARiLexExl01_%s_20161003.mat';
diary('summaryLexExl01.txt')
symbols ={'IBM','NVDA'};
rAll = [];
for iSymbol=1:2
    symbol = symbols{iSymbol};
    inFile = sprintf(fileNameFmt,symbol);
    inFileFull = fullfile(inDir,inFile);
    fprintf('\n%s\nResults for %s\n',repmat('-',1,50),inFileFull)
    titleString=sprintf('%s initial impulse in %s; IRF for %s\n',symbol,shockVar,impactVar);
    load(inFileFull);
    sdw=zeros(0,1);
    is=zeros(0,4);
    resolution = zeros(0,1);
    nmv = length(mva);
    for i=1:nmv
        fprintf('%s %s\n',symbol,repmat('-',1,25))
        mva(i).rwd.rDisplay
        resolution = vertcat(resolution, 1/mva(i).secondFactor);
        sdw = vertcat(sdw,sqrt(mva(i).rwd.var_wAnnual));
        is = vertcat(is,reshape(mva(i).rwd.isGroupBounds',1,[]));
    end
    m=[resolution sdw is];
    colNames = {'resolution' 'sdw' 'min1' 'max1' 'min2' 'max2'};
    disptable(m,colNames,[],'%6.3f')
    
    %   set up table and spreadsheet
    %   symbol, resolution, sdw, min, max, min, max
    vn = {'resolution','sdw','min1','max1','min2','max2'};
    s = string(repmat(symbol,nmv,1));
    rt = horzcat(array2table(s,'VariableNames',{'symbol'}),array2table(m,'VariableNames',vn));
    if isempty(rAll)
        rAll = rt;
    else
        rAll = vertcat(rAll,rt);
    end
    
    [~,iShock]=ismember(shockVar,mva(1).priceNames);
    [~,jVariable]=ismember(impactVar,mva(1).priceNames);
    
    %   Set up graph
    figure('Units','inches','Position',[1,1,6,4]);
    ax = gca;
    ax.XGrid='on';
    ax.XMinorGrid='off';
    ax.FontSize=fontSize;
    ax.TitleFontSizeMultiplier = 1;
    ax.LabelFontSizeMultiplier = 1;
    ax.Title.String=titleString;
    ax.Title.FontWeight='normal';
    ax.YLabel.String='$/share';
    ax.XLabel.String='Time';
    ax.XLimMode = 'manual';
    ax.XGrid='on';
    ax.XMinorGrid='off';
    ax.XScale='log';
    ax.XTick = 10.^(-5:2);
    ax.XLim=[.000001 1300];
    ax.XTickLabel={'10 \mus','100 \mus','1 ms','10 ms','100 ms', '1 sec','10 sec','100 sec'};
    ax.YGrid='on';
    ax.YMinorGrid='off';
    hold
    
    for i=1:length(mva)
        mv = mva(i);
        secondFactor = mv.secondFactor;
        X=mv.irft(2:end)/mv.secondFactor;
        Y=squeeze(mv.irfPacked(jVariable,2:end,iShock));
        semilogx(X,Y);
        switch secondFactor
            case 1; pText='1 sec';
            case 10; pText='100 ms';
            case 100; pText='10 ms';
            case 1000; pText='1 ms';
            case 10000; pText='100 \mus';
            case 100000; pText='10 \mus';
            otherwise; pText=''
        end
        ax = gca;
        dx = 0;
        dy = -(ax.YLim(2)-ax.YLim(1))/30;
        % text(X(1)+dx,Y(1)+dy,[pText],'FontSize',fontSize);
    end
    ytix = arrayfun(@(x) sprintf('$%4.2f',x),ax.YTick,'UniformOutput',false);
    k=ax.YTick<0;
    if any(k); ytix{k}=''; end
    ax.YTickLabel = ytix;
    ax.YLim(1) = min(ax.YLim(1),-0.05);
    ax.YTickMode='manual';  %   'manual' prevents automatic redrawing of ticks if the graph is resized.
%     ax.Children
    h=findall(gca,'Type','Line');
    lc='kbrkbr';
    for i=1:length(h)
        h(i).LineWidth = 1+(length(h)-i)*.2;
        k = mod(i,3);
        switch k
            case 1; h(i).LineStyle='-';
            case 2; h(i).LineStyle=':';
            case 0; h(i).LineStyle='--';
            otherwise; s='-.';
        end
        h(i).Color = lc(i);
        % fprintf('%d %s %f %s\n',i,h(i).LineStyle, h(i).LineWidth, h(i).Color);
    end
    gFile = sprintf('summaryLexExl01 %4s',symbol);
    saveas(gcf,gFile,'fig')
    print('-r300',gFile,'-djpeg')
%     saveas(gcf,gFile,'jpg')
end
tFile = 'summaryLexExl01.xlsx';
writetable(rAll,tFile);
fprintf('summary table written to %s\n',tFile)
diary off