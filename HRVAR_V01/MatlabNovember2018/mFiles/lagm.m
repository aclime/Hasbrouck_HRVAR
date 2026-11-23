function xm = lagm(x, n, f)
%lagm for x a vector returns [x(t) x(t-1) ... x(t-n)]
% f is the fill for initial lagged values (optional, defaults to NaN)
if isvector(x); x=reshape(x,numel(x),1); end
[r,c]=size(x);
if nargin<3; f=NaN; end
xm = NaN(r+n,c*(n+1));
for j=1:c
    for i=0:n
        k=(j-1)*(n+1) + 1;
        xm(:,k+i) = circshift([x(:,j); repmat(f,n,1)],i);
    end
end
xm = xm(1:r,:);
end
