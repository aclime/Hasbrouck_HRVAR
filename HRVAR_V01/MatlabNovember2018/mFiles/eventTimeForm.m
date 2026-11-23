function [xs]=eventTimeForm(sp) 
    %   Return data matrix of prices ordered by time (that is, in event time)
    if ~isa(sp,'spl'); error('Argument should be spl array'); end
    m = length(sp)+1;
    x = NaN(0,m);
    for j=1:length(sp)
        s = NaN(length(sp(j).i),m);
        s(:,1) = sp(j).i;
        s(:,j+1) = sp(j).v;
        x = [x; s];
    end
    n = size(x,1);
    [b,k] = sort(x(:,1));
    xs = NaN(n,m);
    xs = x(k,:);
    for i=2:n
        for j=2:m
            if isnan(xs(i,j)); xs(i,j) = xs(i-1,j); end
        end
    end
    t = xs(:,1);
    j = t(2:end)==t(1:end-1);
    kSameAsNext = [j; false];
    xs = xs(~kSameAsNext,:);
    end