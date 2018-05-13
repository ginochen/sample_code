function var = loadvar_minus1time(varname, Case, factor)
% remove the last time sample from cam cases
var = ncread(sprintf('%s%s.%s.nc',Case.dir,Case.name,varname),varname)*factor;
if ndims(var)==3
   var = var(:,:,end-1);
elseif ndims(var)==4
   var = var(:,:,:,end-1);
else
   error('inconsistent variable dimension')
end
