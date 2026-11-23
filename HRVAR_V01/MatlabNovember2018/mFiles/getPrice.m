function t = getPrice(x,priceVar,timeVar,secondFactor,mDisplay,oneSymbol,oneDate)
%getPrice(inTable,priceVar,timeVar,secondFactor) returns a table: symbol,date,time,price, grouped by secondFactor
% x.Properties
t=[];
if nargin<5; mDisplay=false; end
if nargin>5
    k = strcmp(oneSymbol,x.SYM_ROOT);
    if ~any(k); error(sprintf('ticker symbol %s not found in table',oneSymbol)); end
    x = x(k,:);
end
if nargin>6
    k = int64(oneDate)==int64(x.DATE);
    if ~any(k); error(sprintf('date %d not found in table',oneDate));   end
    x = x(k,:);
end
xDate = x.DATE;
symbol = x.SYM_ROOT;
dates = unique(xDate);
symbols = unique(symbol);
kTime = find(strcmp(x.Properties.VariableNames,timeVar));
if isempty(kTime); error(sprintf('%s not found in table.',timeVar)); end
time = x.(kTime);
kPrice = find(strcmp(x.Properties.VariableNames,priceVar));
if isempty(kPrice); error(sprintf('%s not found in table.',priceVar)); end
price = x.(kPrice);
x2 = table(xDate, symbol, time, price);
x2.Properties.VariableNames{1}='date';
if length(dates)==1 & length(symbols)==1
    x2.time2 = fix(x2.time*secondFactor);
    j = [diff(x2.time2)>0; true]; % last of a run of time stamps.
    t = x2(j,:);
    n1 = height(t);
    % limit to normal hours
    k = t.time>seconds(duration(9,0,0)) & t.time<seconds(duration(17,0,0));
    t = t(k,:);
    %   take only the first time that sets a price ADDED OCTOBER 14, 2018
    k = [true; diff(t.price)~=0];
    t = t(k,:);
    n2 = height(t);
    if mDisplay; fprintf('%d %s %6d observations (%d with nonzero price changes and 9am-5pm)\n',oneDate,oneSymbol,n1,n2); end
else
    nsd = 0;
    for iDate=1:length(dates)
        for iSymbol=1:length(symbols)
            k = dates(iDate)==x2.date & strcmp(symbols(iSymbol),x2.symbol);
            if mDisplay; fprintf('%d %s %6d observations\n',dates(iDate),symbols{iSymbol},sum(k)); end
            x3 = x2(k,:);
            x3.time2 = fix(x3.time*secondFactor);
            j = [diff(x3.time2)>0; true]; % last of a run of time stamps.
            if nsd==0;
                t = x3(j,:);
            else
                t = [t; x3(j,:)];
            end
            nsd = nsd+1;
        end
    end
end
% limit to normal hours
k = t.time>seconds(duration(9,0,0)) & t.time<seconds(duration(17,0,0));
t = t(k,:);
return
end

