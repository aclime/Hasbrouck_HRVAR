classdef xLags < matlab.mixin.Copyable
    %   Defines a range of lags with unrestricted coefficients
    properties
        n = 1   %   lages are x(t-1) ... x(t-n)
        name='xLag'
        kOffset=0;  %   offset within design matrix! NOT 'shift' NOT USED IN xLags ???
        vNames={}
    end
    methods
        function s=toString(o)
            s = sprintf('name=%s; n=%d; kOffset=%d; vNames=',o.name,o.n,o.kOffset);
            s = [s strjoin(o.vNames) '.'];
        end
        function o=polynom(n,name,kOffset)
            if nargin>0; o.n=n; end
            if nargin>1; o.name=name; end
            if nargin>3; o.kOffset = kOffset; end
            for i=1:o.n; o.vNames=[o.vNames, [o.name 'L' int2str(i)]]; end
        end
        function c=designMatrix(o)    %   build design matrix
            c = eye(n);
        end
        function z=sum(o,i1,i2)
            %   returns the sum (across columns) of the design matrix
            %   if i1 and i2 are specified, only columns i1,...,i2 are included
            if nargin==1
                i1 = 1;
                i2 = o.n;
            end
            if i1>i2; error('polynom.sum i1>i2'); end
            z=zeros(1,o.deg+1);
            z(1) = i2-i1+1;
            if o.deg>0; z(2) = (i2+1-i1)*(i2+i1)/2; end
            if o.deg>1; z(3) = (i2+1-i1)*(2*i2*i2 + 2*i1*i1 + 2*i1*i2 + i2 - i1)/6; end
        end
        function ml=maxLength(polya)
            ml = 0;
            for i=1:length(polya)
                ml = max([ml, polya(i).kOffset + polya(i).n]);
            end
        end
        function nt=nTerms(polya)
            nt = 0;
            for i=1:length(polya)
                nt = nt + polya(i).deg+1;
            end
        end
        function su=stateUpdate(polya, state, xIn, xOut)
            % state is nTerms rows x nv cols
            % xIn, xOut are length(polya) rows x nv cols
            su = NaN(size(state));
            k = 0;
            for i=1:length(polya)
                su(k+1,:) = state(k+1,:) + xIn(i,:) - xOut(i,:);
                if polya(i).deg>0; su(k+2,:) = state(k+2,:) + state(k+1,:) + xIn(i,:) - (polya(i).n+1)*xOut(i,:); end
                if polya(i).deg>1; su(k+3,:) = state(k+3,:) + 2*state(k+2,:) + state(k+1,:) + xIn(i,:) - xOut(i,:) * (polya(i).n+1)^2; end
                k = k + polya(i).deg+1;
            end
        end
    end
    methods (Static)
        function c=pSum(x,y,ia,ib,j)
            %   polynomial cross product.
            %   c corresponds to mathematica expression:
            %   Table[Sum[i^d1 (i - j)^d2, {i, ia, ib}], {d1, 0, 2}, {d2, 0, 2}]
            if nargin<2; y=x; end
            if (~isa(x,'polynom')) | (~isa(y,'polynom')); error('polynom.cpSum expects two polynom arguments.'); end
            if nargin<4
                ia = 1;
                ib = min([x.length,y.length]);
            end
            if (ia>ib); error('polynom.cp called with ia>ib.'); end
            if nargin<5; j=0; end
            c = zeros(x.deg+1,y.deg+1);
            for ir=1:x.deg+1
                for ic=1:y.deg+1
                    switch ir
                        case 1
                            switch ic
                                case 1
                                    c(1,1) = 1 + ib - ia;
                                case 2
                                    c(1,2) = -((-1 + ia - ib)*(ia + ib - 2*j))/2;
                                case 3
                                    c(1,3) = -((-1 + ia - ib)*(2*ia^2 + ib + 2*ib^2 + ia*(-1 + 2*ib - 6*j) - 6*ib*j + 6*j^2))/6;
                            end
                        case 2
                            switch ic
                                case 1
                                    c(2,1) = -((-1 + ia - ib)*(ia + ib))/2;
                                case 2
                                    c(2,2) = -((-1 + ia - ib)*(2*ia^2 + ia*(-1 + 2*ib - 3*j) + ib*(1 + 2*ib - 3*j)))/6;
                                case 3
                                    c(2,3) = -((-1 + ia - ib)*(3*ia^3 + ia^2*(-3 + 3*ib - 8*j) + ib*(3*ib^2 + ib*(3 - 8*j) + 2*j*(-2 + 3*j)) + ia*(3*ib^2 - 8*ib*j + 2*j*(2 + 3*j))))/12;
                            end
                        case 3
                            switch ic
                                case 1
                                    c(3,1) = -((-1 + ia - ib)*(-ia + 2*ia^2 + ib + 2*ia*ib + 2*ib^2))/6;
                                case 2
                                    c(3,2) = (-3*ia^4 + ib*(1 + ib)*(3*ib^2 + ib*(3 - 4*j) - 2*j) + 2*ia*j - 3*ia^2*(1 + 2*j) + ia^3*(6 + 4*j))/12;
                                case 3
                                    c(3,3) = (ia - 6*ia^5 - 5*ia*j^2 + 15*ia^4*(1 + j) + 15*ia^2*j*(1 + j) - 10*ia^3*(1 + 3*j + j^2) + ib*(1 + ib)*(-1 + 6*ib^3 + ib^2*(9 - 15*j) + 5*j^2 + ib*(1 - 15*j + 10*j^2)))/30;
                            end
                    end
                end
            end
        end     
    end
end