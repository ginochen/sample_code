function var = loadvar(varname, Case, factor, minusNtime,start,count,stride)
% var dimensions in the ncfile has to be in the order of (time,lev,lat,lon), or
% (time,lat,lon), once read in the order reverses to
% var(lon,lat,lev,time) or var(lon,lat,time)
% 
% start : For an N-dimensional variable, start is a vector of length N of indices specifying the starting location.
% count : Vector of length N specifying the number of elements to read along the corresponding dimensions. 
%         If a particular element of count is Inf, ncread reads data until the end of the corresponding dimension.
% stride : the inter-element spacing along each dimension.
% factor : multiply the var by the factor
% minusNtime : remove the last n time sample from cam cases
% Case.dir : the path for the file 
% Case.name : the name of the file
% varname : the variable name in the file name (must be the same as
% the variable name in the file)
% 
% Example: var(lat,lon,lev,time) has N=4 dimensions, where
% lat=180,lon=360,lev=30,time=3601 so if you want to load last 3
% levels and all the lat's lon's, and subtract the last time index to
% get time=3600 then set 
%
% $ start = [1,1,28,1]; count = [Inf,Inf,3,Inf]; minusNtime=1;
%
if exist('start','var')==0 % no start, count, stride 
   var = ncread(sprintf('%s%s.%s.nc',Case.dir,Case.name,varname),varname)*factor;
elseif exist('start','var') & exist('count','var') & exist('stide','var')==0 % no stride, with start and count only
   var = ncread(sprintf('%s%s.%s.nc',Case.dir,Case.name,varname),varname,start,count)*factor;
elseif exist('start','var') & exist('count','var') & exist('stride','var') % all start, count, stride
   var = ncread(sprintf('%s%s.%s.nc',Case.dir,Case.name,varname),varname,start,count,stride)*factor;
else
   error('input is neither specified, check loadvar code')
end
if ndims(var)==3
   var = var(:,:,1:end-minusNtime);
elseif ndims(var)==4
   var = var(:,:,:,1:end-minusNtime);
elseif ndims(var)==2
else    
   error('inconsistent variable dimension')
end
