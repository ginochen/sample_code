function var = loaddim(varname, Case)
var = ncread(sprintf('%s%s.nc',Case.dir,Case.name),varname);
