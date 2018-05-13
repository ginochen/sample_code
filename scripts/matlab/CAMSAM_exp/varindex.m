
function idv = varindex(str,parm)
  % idv: object variable for variable index
  % example: 
  %   % the index for SST variable
  %   idx.SST = 1
  %   % the str is thus 'SST'
  %   % parm saves the variable in index order
  %   % so suppose two variables 'SST' and 'SLP'
  %   parm = {'SST', 'SLP'}
  for ii=1:numel(str)
    evalc(sprintf('nv = numel(parm.%s)',str{ii}));
    evalc(sprintf('strvar = parm.%s',str{ii}));
    for iv = 1:nv
      evalc(sprintf('idv.%s.%s = %d',str{ii},strvar{iv},iv));
    end
  end
