classdef spd < matlab.mixin.Copyable
    %   spd    sparse vector class suitable for price differences
    %   Before instantiating the first time, use setgetMax to set the maxSize of all the spds in a program.
    %   maxSize is usually related to time. For example, with a sample that could span up to a 24-hour
    %       day, with one observation every minute, maxSize=24 x 60.
    %       With one observation per second, maxSize=24 x 60 x 60, and so on. 
    properties
        i = zeros(1,0)  %indices
        v = zeros(1,0)  %values
        firstValid=[];  %index of first valid observation (usually =i(1)).
        lastValid=[];   %index of last valid observation (defaults to maxSize, but might be the mkt close time.)
        name='spd';     %will usually be set to something more meaningful.
        maxSize=[]      % can be set at time of construction, using the default value from setgetMax
    end
    methods (Static)
        function out=setgetMax(maxSizeIn)
            % setgetMax sets or gets maxSize
            %   example: setget(10^6*3600*24) sets maxSize to the number of microseconds in a day
            persistent maxSize; %   A persistent variable is used in lieu of a Java-type 'static' class variable.
            if nargin; maxSize = maxSizeIn; end
            out = maxSize;
        end
    end
    methods (Access=public)
        function o=spd(i,v,firstValid,lastValid,name)
            %   spd constructor, all arguments are optional
            %   i is a vector of indices
            %   The remaining arguments are optional:
            %   v is a vector of values
            %   1<=firstValid<lastValid<=maxSize
            %   name is a name string or character vector
            if isempty(spd.setgetMax()); error('Attempt to set spd before maxSize is set.'); end
            o.maxSize = spd.setgetMax();
            if nargin==0; return; end
            if nargin>=1
                if ~isvector(i); error('spd(i, v ...): i is not a vector.'); end
                o.i = reshape(i,1,numel(i));
                o.v = ones(1,numel(i));
            end
            if nargin>=2
                if ~isvector(v); error('spd(i, v ...): v (values) is not a vector.'); end
                if length(i)~=length(v); error('spd(i, v ...): i and v not same length.'); end
                o.v = reshape(v,1,numel(i));
                [o.i,p]=sort(o.i);
                o.v=v(p);
            end
            o.firstValid = 1;
            if nargin>=3 && ~isempty(firstValid); o.firstValid=firstValid; end
            o.lastValid = spd.setgetMax();
            if nargin>=4 && ~isempty(lastValid); o.lastValid=lastValid; end
            if nargin==5; o.name=name; end
            o.checkForm
            if any(o.i>o.lastValid)
                fprintf('spd. i has indices > lastValid. (warning)\n')
            end
        end
        function s=toString(o,secondFactor)
            %   toString returns a string with a formatted description.
            %   secondFactor is optional, the number of subperiods in a second
            %   if secondFactor is set, the first and last valid indices are formatted as times.
            s = sprintf('%18s: nnz=%-10d firstValid=%-10d lastValid=%-10d',o.name,length(o.i),o.firstValid,o.lastValid);
            if nargin>1
                prec = log10(secondFactor);
                s = horzcat(s,sprintf(' (%s - %s',timeSec(o.i(1)/secondFactor,prec)),timeSec(o.i(end)/secondFactor,prec),')');
            end
        end
        function sim(s,z)
            %   sim(s,z) fills the vector with random values.
            %   z is either a proportional density (like 0.05) or an absolute number of nonzero elements
            %   the values are random integers between 1 and 100
            if nargin<2; z=0.5; end
            if z>1
                n = z;
            else
                n = floor(z*spd.setgetMax);
            end
            s.i = unique( unidrnd(spd.setgetMax,1,n) );
            s.v = unidrnd(100,size(s.i));
            s.firstValid = 1;
            s.lastValid = spd.setgetMax;
            s.maxSize = spd.setgetMax;
        end
        function checkForm(o)
            %   checkForm performs diagnostic checks
            for k=1:length(o)
                if ~issorted(o(k).i); error('spd i not sorted'); end
                if numel(unique(o(k).i))~=numel(o(k).i); error(['spd ' o(k).name ' has duplicate i''s']); end
                if max(o(k).i)>spd.setgetMax; error('spd i>maxSize'); end
                if ~isempty(find(o(k).i<1, 1)); error('spd i<1'); end
            end
        end
        function rangeCheck(o,first,last,shift)
            % rangeCheck verifies that [first,last] contained in [firstValid+shift,lastValid]
            if nargin==3; shift=0; end
            if shift<0; error('spd.rangeCheck shift<0.'); end
            if first>last; error('spd.rangeCheck. first=%d > last=%d',first,last); end
            if first<o.firstValid+shift; error('spd.rangeCheck. first=%d < firstValid=%d + shift=%d',first,o.firstValid,shift); end
            if last>o.lastValid; error('spd.rangeCheck. last=%d > lastValid=%d',last,o.lastValid); end
        end
        function o2=lag(o,n,xFill)
            %   lag(o,n) returns an spd with the values lagged by n positions;
            if nargin<2; n=1; end
            if nargin<3; xFill=NaN; end
            if n<0; error('spd lag n<0'); end
            o2 = o.copy;
            if n==0; return; end
            if xFill==0
                o2.i = o.i+n;
            else
                o2.i = [1:n o2.i+n];
                o2.v = [NaN(1,n) o2.v];
                o2.firstValid = o2.firstValid+n;
            end
            k = o2.i<=spd.setgetMax();
            o2.i = o2.i(k);
            o2.v = o2.v(k);
            o2.lastValid = min([o2.lastValid+n spd.setgetMax()]);
            o2.name = [o.name '(t-' int2str(n) ')'];
        end
        function n=count(o,first,last,shift)
            %   returns the number of non-zero elements in a given range.
            if nargin<2
                first = o.firstValid;
                last = o.lastValid;
            end
            if nargin<4; shift=0; end
            o.rangeCheck(first,last,shift);
            n = 0.;
            iShift = o.i+shift;
            k = first<=iShift & iShift<=last;
            if isempty(k); return; end
            n = sum(k);
        end
        function s=sum(o,first,last,shift)
            %   sum returns the total of all values in [first,last]
            if nargin<2
                first = o.firstValid;
                last = o.lastValid;
            end
            if nargin<4; shift=0; end
            o.rangeCheck(first,last,shift);
            s=0;
            iShift = o.i+shift;
            k = first<=iShift & iShift<=last;
            if isempty(k); return; end
            s = sum(o.v(k));
        end
        function m=toRow(o)
            %   toRow retuns an spd as a (full) row vector.
            m=NaN;
            if spd.setgetMax>100000; return; end
            m=zeros(1,spd.setgetMax);
            if o.firstValid>1; m(1:o.firstValid-1)=nan; end;
            if o.lastValid<numel(m); m(o.lastValid-1:end)=nan; end;
            m(o.i) = o.v;
        end
        function m=toCol(o)
            %   toCol returns an spd as a (full) column vector.
            m=(o.toRow)';
        end
        function z=dot(x,first,last,xShift,y,yShift)
            %   dot returns the product of two spds over [first,last].
            if ~isa(y,'spd'); error('spd.dot called for non-spd argument.'); end
            z=0;
            if y==x
                k=first<=(x.i+xShift) & (x.i+xShift)<=last;
                if isempty(k); return; end
                z =sum(x.v(k).^2);
                return
            end
            [k,io,iy]=intersect(x.i+xShift,y.i+yShift);
            if isempty(k); return; end
            w=x.v(io) .* y.v(iy);
            w=w(k>=first & k<=last);
            z=sum(w);
        end
    end
end