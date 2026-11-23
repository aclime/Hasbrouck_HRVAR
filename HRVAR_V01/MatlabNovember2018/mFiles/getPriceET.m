function t = getPriceET(x,priceVar,timeVar, inSymbol, inDate)
%getPrice(inTable,priceVar,timeVar,secondFactor) returns a table: symbol,date,time,price, grouped by secondFactor
% x.Properties
t=[];
kTime = find(strcmp(x.Properties.VariableNames,timeVar));
if isempty(kTime); error(sprintf('%s not found in table.',timeVar)); end
time = x.(kTime);
kPrice = find(strcmp(x.Properties.VariableNames,priceVar));
if isempty(kPrice); error(sprintf('%s not found in table.',priceVar)); end
price = x.(kPrice);
date = x.DATE;
symbol = x.SYM_ROOT;
x2 = table(date, symbol, time, price);
k = x2.date==inDate;
x2 = x2(k,:);
k = strcmp(x2.symbol,inSymbol);
x2 = x2(k,:);
k = x2.time>seconds(duration(9,0,0)) & x2.time<seconds(duration(17,0,0));
x2 = x2(k,:);
k = [true; diff(x2.price)~=0];
x2 = x2(k,:);
k = [true; diff(x2.time)~=0];
x2 = x2(k,:);
t = x2;
return
end

