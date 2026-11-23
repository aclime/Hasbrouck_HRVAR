function prices=getPrices(inFile, symbol, v, secondFactor)
%             opts = detectImportOptions(inFile);
%             opts.VariableNames;
x = readtable(inFile);
k = find(strcmp('SYM_ROOT',x.Properties.VariableNames));
x = x(strcmp(symbol,x.(k)),:);
kTimeMid = find(strcmp('timeMid',x.Properties.VariableNames));
tHMS = timeSec(x.(kTimeMid),3);
x = [table(tHMS) x];
kTimeMid = find(strcmp('timeMid',x.Properties.VariableNames));
k = x.(kTimeMid)>[9 45]*[3600; 60] & x.(kTimeMid)<16*3600;  %   9:45-16:00
x = x(k,:);
x.(kTimeMid) = reshape(floor(x.(kTimeMid)*secondFactor),[],1);
if any(diff(x.(kTimeMid))<0); error('time.Mid not sorted in ascending order.'); end
% x.Properties.VariableNames
prices = spl.empty();
for i=1:length(v)
    kv = find(strcmp(v{i},x.Properties.VariableNames));
    k = ~isnan(x.(kv));
    x2 = x(k,:);
    k = [diff(x2.(kTimeMid))>0; true]; % last of a run of time stamps.
    x2 = x2(k,:);
    k = [true; diff(x2.(kv))~=0]; %  first in a run of price changes.
    x2 = x2(k,:);
    p = spl(x2.(kTimeMid),x2.(kv));
    p.name = v{i};
    prices = horzcat(prices, p);
end


