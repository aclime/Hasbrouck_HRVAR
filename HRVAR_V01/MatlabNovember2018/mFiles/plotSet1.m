function plotSet1(fTitle,fontSize)
    ax = gca;
    ax.XGrid='on';
    ax.XMinorGrid='off';
    ax.FontSize=fontSize;
    ax.TitleFontSizeMultiplier = 1;
    ax.LabelFontSizeMultiplier = 1;
    ax.Title.String=fTitle;
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
end

