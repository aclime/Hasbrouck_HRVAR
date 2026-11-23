function timeAxisReset(src,callbackData)
ax = gca;
xVals = 10.^(-6:6);
xLabels = {'1 \mus' '10 \mus' '100 \mus' '1 ms' '10 ms' '100 ms' '1s' '10s' '100s' '1,000s' '10,000s' '100,000s' '1Ms'};
[la,lb]=ismember(ax.XTick, xVals);
ax.XTickLabel=xLabels(lb(la));
end

