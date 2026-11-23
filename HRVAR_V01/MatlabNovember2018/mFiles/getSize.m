function totSize=getSize(this,sDisplay)
if nargin<2; sDisplay=false; end
props = properties(this);
totSize = 0;
for ii=1:length(props)
    currentProperty = getfield(this, char(props(ii)));
    s = whos('currentProperty');
    totSize = totSize + s.bytes;
end
if sDisplay; fprintf('%d bytes\n', totSize); end
end

