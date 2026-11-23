function hmsf = timeSec(t,prec)
%timeString(t,timeFactor,prec): sec past midnight, precision
if nargin==1; prec=0; end
h = floor(t/3600);
m = floor(mod(t,3600)/60);
s = mod(t,60);
s = floor(s*10^prec)*10^(-prec);
if prec==0
    fmt='%02d:%02d:%02d';
else
    fmt = ['%02d:%02d:%0' int2str(prec+3) '.' int2str(prec) 'f'];
end
if isscalar(t)
    hmsf = sprintf(fmt,h,m,s);
else    
    hmsf = arrayfun(@(h,m,s) sprintf(fmt,h,m,s), h,m,s,'UniformOutput',false);
end

