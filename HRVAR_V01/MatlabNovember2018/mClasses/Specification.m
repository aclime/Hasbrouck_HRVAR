classdef Specification 
    %Specification autoregressive or pdl
    enumeration
        AR
        PDL
    end
    methods
        function o=Specification()
        end
        function r=isAR(o)
            r = o==Specification.AR;
        end
        function r=isPDL(o)
            r = o==Specification.PDL;
        end
    end
end

