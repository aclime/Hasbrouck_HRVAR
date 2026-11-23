classdef polynom < matlab.mixin.Copyable
    %   The (implicit) object here is a design matrix of deg+1 rows and n columns
    %   The rows correspond to degree=0,1,...,deg polynomial degrees.
    %   THE CODE HERE IS ONLY IMPLEMENTED for deg<=2
    %   Many of the routines are implemented for polynom arrays (identified by first argument: "polya")
    properties
        deg=0       %   degree of polynomial
        n=1         %   length of polynomial
        name='poly'
        kOffset=0   %   offset within design matrix! NOT 'shift'
        vNames={}   %   coefficient names (cell array of character vectors)
    end
    methods
        function s=toString(o)
            s = sprintf('name=%s; deg=%d; n=%d; kOffset=%d; vNames=',o.name,o.deg,o.n,o.kOffset);
            s = [s strjoin(o.vNames) '.'];
        end
        function o=polynom(deg,n,name,kOffset)
            if nargin>0; o.deg=deg; end
            if nargin>1; o.n=n; end
            if nargin>2 && ~isempty(name); o.name=name; end
            if nargin>3; o.kOffset = kOffset; end
            if o.deg>2; error('polynom degrees>2 aren''t supported'); end
            for i=0:o.deg; o.vNames=[o.vNames, [o.name 'd' int2str(i)]]; end
        end
        function vNamesSet(o)   %   used to recompute vNames after the degree has changed.
            if o.deg>2; error('polynom degrees>2 aren''t supported'); end
            o.vNames={};
            for i=0:o.deg; o.vNames=[o.vNames, [o.name 'd' int2str(i)]]; end
        end
        function c=designMatrix(o)    %   build design matrix
            c = zeros(o.nTerms, o.maxLength);
            kr = 0;
            for i=1:length(o);
                krEnd = kr + o(i).deg+1;
                c(kr+1:krEnd, o(i).kOffset+(1:o(i).n)) = repmat(1:o(i).n,o(i).deg+1,1).^transpose(0:o(i).deg);
                kr = krEnd;
            end
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
        function ml=maxLength(polya)    %   maximum lag of a set (array) of polynomials
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
            %   where d1 and d2 correspond to the max degrees of x and y.
            if nargin<2; y=x; end
            if nargin<4
                ia = 1;
                ib = min([x.length,y.length]);
            end
            if nargin<5; j=0; end
            if (~isa(x,'polynom')) | (~isa(y,'polynom')); error('polynom.pSum expects two polynom arguments.'); end
            if (ia>ib); error('polynom.cp called with ia>ib.'); end
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