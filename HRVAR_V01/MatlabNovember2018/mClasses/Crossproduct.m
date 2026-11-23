classdef Crossproduct
    % Crossproduct is a container for sparse crossproduct routines.
    % Polynoms (including arrays), spd's, and spl's
    % Modification for polynom arrays April 24 2017 (used in MVARd).
    % Modification for non-polynoms January, 2018.
    % In these routines, xShift and yShift are the lags relative to the positions given in the data array.
    methods (Static)
        %%  crossproducts of lagm and ... (used in AR estimations)
        %Conceptually: lagm for x an spd vector is [x(t) x(t-1) ... x(t-n)] with NaNs in missing places
        function z=lagmxconst(x,xLagMax,first,last)
            if ~isa(x,'spd'); error('Crossproduct.lagmxconst(x,xLagMax,first,last) x is not an spd object.'); end
            z = zeros(xLagMax+1,1);
            for iShift=0:xLagMax
                z(iShift+1) = Crossproduct.spdxconst(x,iShift,first,last);
            end
        end
        function z=lagmxspd(x,xLagMax,y,yShift,first,last)
            if ~isa(x,'spd') | ~isa(y,'spd'); error('Crossproduct.lagmxspd(x,xLagMax,y,yShift,first,last) x or y not spd object.'); end
            z = zeros(xLagMax+1,1);
            for iShift=0:xLagMax+1
                z(iShift+1) = Crossproduct.spdxspd(x,iShift,y,yShift,first,last);
            end
        end
        function z=lagmxspl(x,xLagMax,y,yShift,first,last)
            if ~isa(x,'spd'); error('Crossproduct.lagmxspl(x,xLagMax,y,yShift,first,last) x is not an spd object.'); end
            if ~isa(y,'spl'); error('Crossproduct.lagmxspl(x,xLagMax,y,yShift,first,last) y is not an spl object.'); end
            z = zeros(xLagMax+1,1);
            for iShift=0:xLagMax
%                 z(iShift+1) = Crossproduct.spdxspl(x,iShift,y,yShift,first,last);
                z(iShift+1) = Crossproduct.splxspd(y,yShift,x,iShift,first,last);
            end
        end
        function z=lagmxlagm(x,xLagMax,y,yLagMax,first,last)
            if ~isa(x,'spd'); error('Crossproduct.lagmxlagm(x,xLagMax,y,yLagMax,first,last) x is not an spd object.'); end
            if ~isa(y,'spd'); error('Crossproduct.lagmxlagm(x,xLagMax,y,yLagMax,first,last) y is not an spd object.'); end
            z = zeros(xLagMax+1,yLagMax+1);
            for iShift=0:xLagMax
                for jShift=0:yLagMax
                    z(iShift+1,jShift+1) = Crossproduct.spdxspd(x,iShift,y,jShift,first,last);
                end
            end
        end
        %%  crossproducts of polynom and ...    
        function z=polyxpoly2(x,px,xShift,y,py,yShift,first,last) % SET UP FOR ARRAYS
            %   returns an x.deg+1 x y.deg+1 matrix
            %   works with px and py arrays
            if nargin~=8; error('Crossproduct.polyxpoly expects eight arguments.'); end
            if ~isa(px,'polynom') | ~isa(py,'polynom'); error('Crossproduct.polyxpoly expects two polynom arguments.'); end
            if ~isa(x,'spd') | ~isa(y,'spd'); error('Crossproduct.polyxpoly expecting two spd arguments.'); end
            if xShift<0 | yShift<0; error('Crossproduct.polyxpoly called with negative xShift/yShift.'); end
            nx = 0;
            for i=1:length(px); nx = nx + px(i).deg + 1; end
            ny = 0;
            for i=1:length(py); ny = ny + py(i).deg + 1; end
            z = zeros(nx,ny);
            kpx = 0;
            for ip=1:length(px)
                kpxEnd = kpx + px(ip).deg + 1;
                xi = x.i + xShift + px(ip).kOffset;
                kpy = 0;
                for jp=1:length(py)
                    yi = y.i + yShift + py(jp).kOffset;
                    zz = zeros(px(ip).deg+1, py(jp).deg+1);
                    kpyEnd = kpy + py(jp).deg + 1;
                    ky0Last = 0;
                    ky0 = 1;
                    ky0 = find(yi+py(jp).n-1>=first,1);
                    kyLast = find(yi<=last,1,'last');
                    nPrint = 0;
                    kxFirst = find(xi+px(ip).n-1>=first,1);
                    kxLast = find(xi<=last,1,'last');
%                     for kx=1:length(xi)
                    for kx=kxFirst:kxLast
%                         if xi(kx)+px(ip).n-1<first; continue; end
%                         if xi(kx)>last; break; end
                        firstNZ = false;
%                         if (kx==kxFirst || kx==kxLast || (ky0>ky0Last && nPrint<=100 && mod(kx,10)==0))
%                             fprintf('ip=%d jp=%d ky0=%d ky0Last=%d kx=%d\n',ip,jp,ky0,ky0Last,kx); 
%                             nPrint = nPrint+1;
%                         end
                        ky0Last = ky0;
%                         for ky=ky0:length(yi)
                        for ky=ky0:kyLast
%                             if yi(ky)+py(jp).n-1<first
% %                                 ky0 = ky;
%                                 continue; 
%                             end
%                             if yi(ky)>last | yi(ky)>xi(kx)+px(ip).n-1; break; end;
                            if yi(ky)>xi(kx)+px(ip).n-1; break; end
                            pcp =   Crossproduct.crossSumPolyxPoly(xi(kx), yi(ky), px(ip), py(jp), first, last);
%                             pcpz = Crossproduct.crossSumPolyxPoly0(xi(kx), yi(ky), px(ip), py(jp), first, last);
%                             if any(pcp~=pcpz)
%                                 pcp
%                                 pcpz
%                                 error(sprintf('pcp~=pcpz: %15.5e\n',pcp,pcpz)); 
%                             end
                            zz = zz + x.v(kx)*y.v(ky)*pcp;
%                         if ~firstNZ && any(pcp~=0,'all')  %   This works on local machine running R2018b, but not HPC using R2018a
                            if ~firstNZ && any(any(pcp~=0))
%                             if ~firstNZ && pcpz~=0
                                firstNZ = true;
                                ky0 = ky;   %   starting point for next pass through outer loop.
                            end
                        end
                    end
                    z(kpx+1:kpxEnd, kpy+1:kpyEnd) = zz;
                    kpy = kpyEnd;
                end
                kpx = kpxEnd;
            end
        end
        function z=polyxpoly(x,px,xShift,y,py,yShift,first,last,displayLevel) % SET UP FOR ARRAYS
            %   returns an x.deg+1 x y.deg+1 matrix
            %   works with px and py arrays
            if nargin<9; displayLevel=0; end
            if nargin<8; error('Crossproduct.polyxpoly expects at least eight arguments.'); end
            if ~isa(px,'polynom') | ~isa(py,'polynom'); error('Crossproduct.polyxpoly expects two polynom arguments.'); end
            if ~isa(x,'spd') | ~isa(y,'spd'); error('Crossproduct.polyxpoly expecting two spd arguments.'); end
            if xShift<0 | yShift<0; error('Crossproduct.polyxpoly called with negative xShift/yShift.'); end
            if displayLevel>0; fprintf('polyxpoly called for x=%s (shift=%d) and y=%s (shift=%d)\n',x.name,xShift,y.name,yShift); end
            cTic=tic;
            z=[];
            for ip=1:length(px)
                if displayLevel>1; fprintf('px: %s\n',px(ip).toString); end
                xi = x.i + xShift + px(ip).kOffset;
                kxFirst = find(xi+px(ip).n-1>=first,1);
                kxLast = find(xi<=last,1,'last');
                zr = [];
                for jp=1:length(py)
                    if displayLevel>1; fprintf('py: %s\n',py(jp).toString); end
                    yi = y.i + yShift + py(jp).kOffset;
                    zz = zeros(px(ip).deg+1, py(jp).deg+1);
                    if displayLevel>2; tic; end;
                    for kx=kxFirst:kxLast
                        kyFirst = find(yi+py(jp).n-1>=xi(kx),1);
                        kyLast = find(yi<=xi(kx)+px(ip).n-1,1,'last');
                        for ky=kyFirst:kyLast
                            pcp =   Crossproduct.crossSumPolyxPoly(xi(kx), yi(ky), px(ip), py(jp), first, last);
                            zz = zz + x.v(kx)*y.v(ky)*pcp;
                        end
                    end
                    if displayLevel>2; toc; end;
                    zr = horzcat(zr,zz);
                end
                z = vertcat(z,zr);
            end
            if displayLevel>0; toc(cTic); end
        end
        function z=polyxspl(x,px,xShift,y,yShift,first,last) %  ARRAYS
            %   returns a px.deg+1 column vector
            if nargin~=7; error('Crossproduct.polyxspl expects seven arguments.'); end
            if ~isa(px,'polynom'); error('Crossproduct.polyxspl expects one polynom argument.'); end
            if ~isa(x,'spd') | ~isa(y,'spl'); error('Crossproduct.polyxspl expects one spd and one spl argument.'); end
            if xShift<0 | yShift<0; error('Crossproduct.polyxspl called with negative xShift/yShift.'); end
            yi = [y.i inf] + yShift;
            yiEnd = yi(2:end)-1;
            yi2 = min([yi; repmat(last+1,1,length(yi))],[],1);
            yn = diff(yi2);
            z = [];
            py = polynom(0,1,'py'); %   This is an artificial device to propogate the price level.
            for ip=1:length(px)
                zz = zeros(px(ip).deg+1,1);
                xi = x.i + xShift + px(ip).kOffset;
                kxFirst = find(xi+px(ip).n-1>=first,1);
                kxLast = find(xi<=last,1,'last');
                for kx=kxFirst:kxLast
                    kyFirst = find(yiEnd>=xi(kx),1);
                    kyLast = find(yi<=xi(kx)+px(ip).n-1,1,'last');
                    for ky=kyFirst:kyLast
                        py.n = yn(ky);
%                         m = min([yi(ky+1),last+1]) - yi(ky);
%                         py.n = m;
%                         if yi(ky)+py.n-1<first; continue; end
%                         if yi(ky)+py.n-1<xi(kx); continue; end; % new?
%                         if yi(ky)>last | yi(ky)>xi(kx)+px(ip).n-1; break; end;
                        pcp = Crossproduct.crossSumPolyxPoly(xi(kx), yi(ky), px(ip), py, first, last);
                        zz = zz + x.v(kx)*y.v(ky)*pcp;
                    end
                end
                z = vertcat(z,zz);
            end
        end
        function z=polyxspd(x,px,xShift,y,yShift,first,last)
            %   returns a column vector with px.deg+1 rows
            if nargin~=7; error('Crosspropoduct.polyxspd expects seven arguments.'); end
            if ~isa(px,'polynom'); error('Crossproduct.polyxspd expects one polynom argument.'); end
            if ~isa(x,'spd') | ~isa(y,'spd'); error('Crossproduct.polyxspd expecting two spd arguments.'); end
            if xShift<0 | yShift<0; error('Crossproduct.polyxspd called with negative xShift/yShift.'); end
            yi = y.i + yShift;
            z = zeros(0,1);
            for ip=1:length(px) %   Loop over px array
                zz = zeros(px(ip).deg + 1, 1);
                xi = x.i + xShift + px(ip).kOffset;
                kxFirst = find(xi+px(ip).n-1>=first,1);
                kxLast = find(xi<=last,1,'last');
                for kx=kxFirst:kxLast
                    kyFirst = find(yi>=xi(kx),1);
                    kyLast = find(yi<=xi(kx)+px(ip).n-1,1,'last');
                    for ky=kyFirst:kyLast
                        pcp = Crossproduct.crossSumPolyxspd(xi(kx), yi(ky), px(ip), first, last);
                        zz = zz + x.v(kx)*y.v(ky)*pcp';
                    end
                end
                z = vertcat(z,zz);
            end
        end
        function z=polyxspdOld(x,px,xShift,y,yShift,first,last)
            %   returns a column vector with px.deg+1 rows
            if nargin~=7; error('Crosspropoduct.polyxspd expects seven arguments.'); end
            if ~isa(px,'polynom'); error('Crossproduct.polyxspd expects one polynom argument.'); end
            if ~isa(x,'spd') | ~isa(y,'spd'); error('Crossproduct.polyxspd expecting two spd arguments.'); end
            if xShift<0 | yShift<0; error('Crossproduct.polyxspd called with negative xShift/yShift.'); end
            nxPrintMax=3;
            nyPrintMax=100;
            nr = sum(arrayfun(@(x) x.deg+1, px))
            yi = y.i + yShift;
            z = zeros(nr,1);
            kx0 = 0;
            for ip=1:length(px)
                kxEnd = kx0 + px(ip).deg + 1;
                zz = zeros(px(ip).deg + 1, 1);
                xi = x.i + xShift + px(ip).kOffset;
                ky0 = find(yi>=first,1);
                ky0Last = 0;
                kyLast = find(yi<=last,1,'last');
                kxFirst = find(xi+px(ip).n-1>=first,1);
                kxLast = find(xi<=last,1,'last');
                nxPrint = 0;
                for kx=kxFirst:kxLast
%                     if xi(kx)+px(i).n-1<first; continue; end
%                     if xi(kx)>last; break; end
                    if nxPrint<=nxPrintMax; fprintf('ip=%d kx=%d xi(kx)=%d\n', ip, kx, xi(kx)); end
                    nxPrint =nxPrint+1;
                    if nxPrint>nxPrintMax; return; end
                    firstNZ = false;
%                     if (kx==kxFirst || kx==kxLast || (ky0>ky0Last && nPrint<=100 && mod(kx,10)==0))
%                         fprintf('ip=%d ky0=%d ky0Last=%d kx=%d\n',ip,ky0,ky0Last,kx);
%                         nPrint = nPrint+1;
%                     end
                    nyPrint = 0;
                    ky0Last = ky0;
                    
                    for ky=ky0:kyLast
                        if yi(ky)<first-100; continue; end
                        if nyPrint<=nyPrintMax; fprintf('nyPrint=%d ky0=%d ky=%d yi(ky)=%d\n',nyPrint,ky0,ky,yi(ky)); end
                        nyPrint = nyPrint+1;
                        
                        %                         if yi(ky)<first; continue; end
                        %                         if yi(ky)>last
%                         s=sprintf('kx=%d xi(kx)=%d ky=%d yi(ky)=%d ky0=%d\n',kx,xi(kx),ky,yi(ky),ky0);
                        if yi(ky)>xi(kx)+px(ip).n-1; break; end
                        pcp = Crossproduct.crossSumPolyxspd(xi(kx), yi(ky), px(ip), first, last);
%                         if any(pcp~=0); fprintf('     ky0=%d ky=%d yi(ky)=%d\n',ky0,ky,yi(ky)); disp(pcp); end
%                         fprintf('ky=%d pcp=%d\n',ky,pcp);
                        zz = zz + x.v(kx)*y.v(ky)*pcp';
%                         if ~firstNZ && any(pcp~=0,'all')  %   This works on local machine running R2018b, but not HPC using R2018a
                        if ~firstNZ && any(any(pcp~=0))
                            firstNZ = true;
                            ky0 = ky;   %   starting point for next pass through outer loop.
                        end
                    end
                end
                z(kx0+1:kxEnd,1) = zz;
                kx0 = kxEnd;
            end
        end
        
        function z=polyxconst(x,px,xShift,first,last)   % ARRAYS
            %   returns a row vector of length px.deg+1
            %   px can be an array
            nc = 0;
            for i=1:length(px); nc=nc+px(i).deg+1; end
            z = zeros(1,nc);
            kz = 0;
            for i=1:length(px)
                kzEnd = kz+px(i).deg++1;
                zz = zeros(1,px(i).deg+1);
                xi = x.i+xShift+px(i).kOffset;
                for kx=1:length(xi)
                    if xi(kx)+px(i).n-1<first; continue; end
                    if xi(kx)>last; break; end
                    ka = max([1,first-xi(kx)+1]);
                    kb = min([px(i).n,last-xi(kx)+1]);
                    if ka>kb; continue; end
                    zz = zz + x.v(kx)*sum(px(i),ka,kb);
                end
                z(kz+1:kzEnd) = zz;
                kz = kzEnd;
            end
        end
        %%  crossproducts of spl and ...
        function z=splxspl(x,xShift,y,yShift,first,last)
            z = 0;
            if ~isa(x,'spl') | ~isa(y,'spl'); error('Crossproduct.splxspl expects two spl arguments.'); end
            [iBoth,i1,i2] = Crossproduct.merge(x,xShift,y,yShift);
            iBoth = [iBoth inf];
            if iBoth(1)>first; 
                fprintf('x: %12s i=%12d %12d ...\n',x.name,x.i(1),x.i(2));
                fprintf('y: %12s i=%12d %12d ...\n',y.name,y.i(1),y.i(2));
                error('Crossproduct.splxspl iBoth(1)=%d > first=%d.',iBoth(1),first); 
            end
            for i=1:length(iBoth)-1
                if iBoth(i)>last; break; end
                if iBoth(i+1)<=first; continue; end
                k1 = max(iBoth(i),first);
                k2 = min(iBoth(i+1)-1,last);
                d = k2-k1+1;
                j1 = i1(i);
                j2 = i2(i);
                z = z + x.v(j1)*y.v(j2)*d;
            end
        end
        function z=splxspd(x,xShift,y,yShift,first,last)
            z = 0;
            if ~isa(x,'spl'); error('Crossproduct.splxspd expects one spl argument.'); end
            if ~isa(y,'spd'); error('Crossproduct.splxspd expects one spd argument.'); end
            xi = x.i+xShift;
            yi = y.i+yShift;
            for i=1:length(yi)
                t = yi(i);
                if t>last; break; end
                if t<first; continue; end
                k = find(xi<=t,1,'last');
                if ~isempty(k); z = z + y.v(i)*x.v(k); end;
            end
        end
        function z=splxconst(x,xShift,first,last)
            if ~isa(x,'spl'); error('Crossproduct.splxconst expecting an spl argument.'); end
            z = 0;
            xi = [x.i inf]+xShift; %   append to o.i to simplify looping
            for k=1:length(xi)
                if xi(k)>last; return; end
                if xi(k+1)<first; continue; end
                w = min([xi(k+1),last+1]) - max([xi(k),first]);
                z = z + x.v(k)*w;
            end
        end
        %%  crossproducts of spd and ...
        function z=spdxspd(x,xShift,y,yShift,first,last)
            z=0;
            xi = x.i+xShift;
            yi = y.i+yShift;
            [k,ix,iy]=intersect(xi,yi);
            if isempty(k); return; end
            w=x.v(ix) .* y.v(iy);
            w=w(k>=first & k<=last);
            z=sum(w);
        end
        function z=spdxconst(x,xShift,first,last)
            if nargin~=4; error('Crossproduct.spdxconst expecting four arguments.'); end
            z = 0;
            xi = x.i + xShift;
            k = xi>=first & xi<=last;
            if ~isempty(k)
                z = sum(x.v(k));
            end
        end
%%  Helper routines
        function [iResult,ix,iy]=merge(x,xShift,y,yShift)
            xi = x.i + xShift;
            yi = y.i + yShift;
            iBoth = union(xi,yi);
            [~,ix]=ismember(iBoth,xi);
            [~,iy]=ismember(iBoth,yi);
            if ix(1)~=0 & iy(1)~=0
                iFirst=1;
            else
                iFirst=0;
            end
            for i=2:length(iBoth);
                if ix(i)==0; ix(i)=ix(i-1); end
                if iy(i)==0; iy(i)=iy(i-1); end
                if iFirst==0 & ix(i)~=0 & iy(i)~=0; iFirst=i; end
            end
            iResult = iBoth(iFirst:end);
            ix = ix(iFirst:end);
            iy = iy(iFirst:end);
        end
        function c=crossSumPolyxPoly(i,j,x,y,firstValid,lastValid)  %   ARRAYS
            %   For the first series/polynomial:
            %   i = ii + a where ii is the position of the point in the original data (not 'lagged')
            %       and a is the lag offset.
            %   j = jj + b ... similarly
            %   x and y are the associated polynoms.
            %   firstValid and lastValid specify the valid range for computation.
            %   routine is set up for polynomial arrays, but does not consider offsets.
            %   (all polynomials are presumed to start at the same place)
            if ~isa(x,'polynom') || ~isa(y,'polynom'); error('Crossproduct.crossSumPolyxPoly expects x and y to be polynoms.'); end
            nx = 0;
            for ip=1:length(x); nx = nx + x(ip).deg+1; end
            ny = 0;
            for ip=1:length(y); ny = ny + y(ip).deg+1; end
            c = zeros(nx,ny);
            kx =0;
            for ip=1:length(x);
                kxEnd = kx + x(ip).deg+1;
                ky=0;
                for jp=1:length(y)
                    kyEnd = ky+y(jp).deg+1;
                    k = j-i;
                    ka = max([1, k+1, firstValid-i+1]);
                    kb = min([x(ip).n, y(ip).n+k, lastValid-j+k+1]);
                    if ka>kb; continue; end
                    cij = polynom.pSum(x(ip),y(jp),ka,kb,k);
                    c(kx+1:kxEnd, ky+1:kyEnd) = cij;
                    ky = kyEnd;
                end
                kx = kxEnd;
            end  
        end
        function c=crossSumPolyxPoly0(i,j,x,y,firstValid,lastValid)
            %   degree zero. no argument checking.
            c = 0;
            k = j-i;
            ka = max([1, k+1, firstValid-i+1]);
            kb = min([x.n, y.n+k, lastValid-j+k+1]);
            if ka>kb; return; end
            c = polynom.pSum(x,y,ka,kb,k);
        end
        function c=crossSumPolyxspd(i,j,x,firstValid,lastValid)
            %   cross-product polynom vs. one-position spd
            if ~isa(x,'polynom'); error('Crossproduct.crossSumPolyxspd expects x to be a polynom.'); end
            c = zeros(1,x.deg+1);
            k = j-i;
            ka = max([1, k+1, firstValid-i+1]);
            kb = min([x.n, k+1, lastValid-j+k+1]);
            if ka>kb; return; end
            c = x.sum(ka,kb);
        end
    end
end

