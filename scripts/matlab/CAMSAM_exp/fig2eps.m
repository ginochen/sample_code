function fig2eps(path,fname) 
if ~exist('fname') % having included single file in this function, fix it later when necessary
  tmp = dir([path '/*.fig']);
  fname = {tmp.name};
end
for ifi = 1:numel(fname)
%  set(0, 'CurrentFigure', 1);
%  clf reset
  [~, tmp] = fileparts(fname{ifi})
  eval(sprintf('open %s.fig',tmp))
  eval(sprintf('print -depsc %s',tmp))
end
%eval(sprintf('print -dpng %s',fname))
