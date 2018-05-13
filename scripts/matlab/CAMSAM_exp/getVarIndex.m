function idv = getVarIndex(str,parm)
% Purpose: get the variable index of varCAM0D varCAM1D, ..., etc from parm
% str: variable name
% idv: variable index
for ii=1:numel(str)
   evalc(sprintf('nv = numel(parm.%s)',str{ii}));
   evalc(sprintf('strvar = parm.%s',str{ii}));
   for iv = 1:nv
      evalc(sprintf('idv.%s.%s = %d',str{ii},strvar{iv},iv));
   end
end
