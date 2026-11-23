function d = etDuration(et)
d = seconds(et);
s = seconds(d);
% d.Format='mm:ss.SS';
% return
% fprintf('%d\n',d)
if s<10
    d.Format='s';
elseif s<600
    d.Format='mm:ss.SS';
elseif s<3600
    d.Format='mm:ss';
else
    d.Format='hh:mm:ss';
end
end
