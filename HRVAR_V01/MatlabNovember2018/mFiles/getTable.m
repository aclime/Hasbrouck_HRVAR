function x=getTable(inFile, symbol, secondFactor)
%             opts = detectImportOptions(inFile);
%             opts.VariableNames;
x = readtable(inFile);
x = x(strcmp(symbol,x.SYM_ROOT),:);
tHMS = timeSec(x.timeMid,3);
x = [table(tHMS) x];
k = x.timeMid>[9 45]*[3600; 60] & x.timeMid<16*3600;  %   9:45-16:00
x = x(k,:);
x.timeMid = reshape(floor(x.timeMid*secondFactor),[],1);
if any(diff(x.timeMid)<0); error('timeMid not sorted in ascending order.'); end

