function disptableLog(x,rowNames,n,xBase)
%disptableLog(x,rowNames,n,xBase)
%   displays starting, ending, and representative intermediate columns of x
%   xBase is the starting column index (usually 0 or 1. usually 0 for irfs)
if nargin<4; xBase=1; end
if nargin<3; n=5; end
x = squeeze(x);
if isvector(x); x=reshape(x,1,[]); end
ll = size(x,2);
m=ceil(log10(ll));
k = unique( horzcat(reshape( (xBase:n)'*(10.^(0:m)), 1, []), ll+xBase-1) );
k = k(k<=ll+xBase-1);
s = strsplit(int2str(k));
j = k+1-xBase;
disptable(x(:,j),s,rowNames);
end

