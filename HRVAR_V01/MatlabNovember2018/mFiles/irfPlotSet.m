function irfPlotSet(irfTitle,fontSize)
ax = gca;
ax.XGrid='on';
ax.XMinorGrid='off';
ax.FontSize=fontSize;
ax.Title.String=irfTitle;
ax.Title.FontWeight='normal';
ax.YLabel.String='$/share';
ax.XLabel.String='Time';
% ax.XLimMode = 'manual';
ax.XGrid='on';
ax.XMinorGrid='off';
ax.Title.FontWeight='normal';
ax.YGrid='on';
ax.YMinorGrid='off';
end


