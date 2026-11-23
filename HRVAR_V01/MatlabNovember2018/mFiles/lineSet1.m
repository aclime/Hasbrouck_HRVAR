function lineSet1
h=findall(gca,'Type','Line');
lc='kbrkbr';
for i=1:length(h)
    h(i).LineWidth = 1+(length(h)-i)*.2;
    k = mod(i,3);
    switch k
        case 1; h(i).LineStyle='-';
        case 2; h(i).LineStyle=':';
        case 0; h(i).LineStyle='--';
        otherwise; h(i).LineStyle='-.';
    end
    h(i).Color = lc(i);
    % fprintf('%d %s %f %s\n',i,h(i).LineStyle, h(i).LineWidth, h(i).Color);
end
end

