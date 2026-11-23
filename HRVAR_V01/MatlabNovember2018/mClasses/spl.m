classdef spl < spd
    %   spl is a sparse matrix representing persistent levels.
    %   In a regular sparse matrix, empty elements are assumed to be zero.
    %   This works for price differences, for not for price levels (or
    %   price differences). For example the difference between bid1 and
    %   bid2 (in alternate markets) might be 1 (persisting for five
    %   minutes) or zero persisting for five minutes. In the spl class, the
    %   price is presumed to persist until it is explicitly changed.
    %   The spl class adds no new properties to spd.
    %   The difference between spd and spl is one of interpretation.
    methods (Access=public)
        function o=spl(i,v,firstValid,lastValid,name)
            %   spl constructor
            %   i and v are indices and values
            %   firstValid, lastValid and name are optional
            o = o@spd();
            o.name = 'spl';
            if nargin==0; return; end
            if nargin>=1
                if ~isvector(i); error('spl(i, v ...): i is not a vector.'); end
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
        end
        function z=toCol(o)
            %   toCol returns a full column vector
            z = NaN(o.lastValid,1);
            z(o.i) = o.v;
            for i=o.firstValid:o.lastValid
                if isnan(z(i)) & i>=2; z(i)=z(i-1); end
            end
        end
        function z=toRow(o)
            %   toRow returns a full row vector
            z = toCol(o)';
        end
        function s=sum(o,first,last,shift)
            if nargin<2
                first = o.firstValid;
                last = o.lastValid;
            end
            if nargin<4; shift = 0; end
            s = 0;
            noi = length(o.i);
            oi = [o.i inf]+shift; %   append to o.i to simplify looping
            for k=1:noi
                if oi(k)>last; return; end
                if oi(k+1)<first; continue; end
                w = min([oi(k+1),last+1]) - max([oi(k),first]);
                s = s + o.v(k)*w;
            end
        end
        function s=diff(o)
            o.checkForm
            v = diff(o.v);
%             s = spd;
%             s.set(o.i(2:end),v,o.firstValid+1,o.lastValid,['d' o.name]);
            s = spd(o.i(2:end), v, o.firstValid+1, o.lastValid, ['d' o.name]);
        end
        function o=splMinus(o1,o2)
            [iResult,i1,i2]=merge(o1,o2);
            v = o1.v(i1)-o2.v(i2);
            o = spl;
            o.i = iResult;
            o.v = v;
            if ~isempty(o1.firstValid) && ~isempty(o2.firstValid)
                o.firstValid=max([o1.firstValid o2.firstValid]);
            end
            if ~isempty(o1.lastValid) && ~isempty(o2.lastValid)
                o.lastValid = min([o1.lastValid o2.lastValid]);
            end
            o.name = [o1.name '-' o2.name];
        end
        function o=splTimes(o1,o2)
            [iResult,i1,i2]=merge(o1,o2);
            v = o1.v(i1)*o2.v(i2);
            o = spl(iResult,v);
        end
        function [iResult,i1,i2]=merge(o1,o2)
            %   merge is used to locate indices where one or both of the spls changes.
            iBoth = union(o1.i,o2.i);
            [~,i1]=ismember(iBoth,o1.i);
            [~,i2]=ismember(iBoth,o2.i);
            iFirst=0;
            if i1(1)~=0 & i2(1)~=0; iFirst=1; end
            for i=2:length(iBoth);
                if i1(i)==0; i1(i)=i1(i-1); end
                if i2(i)==0; i2(i)=i2(i-1); end
                if iFirst==0 & i1(i)~=0 & i2(i)~=0; iFirst=i; end
            end
            iResult = iBoth(iFirst:end);
            i1 = i1(iFirst:end);
            i2 = i2(iFirst:end);
            %   Check
            if any(i1==0); 
                k = find(i1==0,1);
                fprintf('spl.merge i1 zero check at position %d. i2(.)=\n',k,i2(k)); 
            end
            if any(i2==0); 
                k = find(i2==0,1);
                fprintf('spl.merge i2 zero check at position %d. i1(.)=\n',k,i1(k)); 
            end
        end
        function z=dotSpl(s1,first,last,s2)
            %dotSpl(s1,first,last,s2) returns the dot product of s1 and s2
            %
            z = 0;
            if ~isa(s2,'spl'); error('spl. dotSpl called with a non-spl argument'); end
            [iBoth,i1,i2] = s1.merge(s2);
            iBoth = [iBoth inf];
            for i=1:length(iBoth)-1
                if iBoth(i)>last; break; end
                if iBoth(i+1)<=first; continue; end
                k1 = max(iBoth(i),first);
                k2 = min(iBoth(i+1)-1,last);
                d = k2-k1+1;
                j1 = i1(i);
                j2 = i2(i);
                z = z + s1.v(j1)*s2.v(j2)*d;
            end
        end
        function z=spdDot(o,s,firstValid,lastValid)
            z = 0;
            if ~isa(s,'spd'); error('spl: spdDot not called with an spd object.'); end
            for i=1:length(s.i)
                t = s.i(i);
                if t<firstValid | t>lastValid; continue; end;
                k = find(o.i<=t,1,'last');
                if ~isempty(k); z = z + s.v(i)*o.v(k); end;
            end
        end
    end
end


