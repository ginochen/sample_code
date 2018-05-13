function iivar = var_index(var)
if ndims(var) ~= 3
   error('dimension has to be 3, (lon,lat,time)')
end
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/index/basin_index.mat'); 
% tropics non-zero PRECCDZM index 20S < lat < 20N, basin-wise
dim = size(var);
ngp = dim(1)*dim(2); %number of grid points
for ib = 1:numel(basin)
   iivar{ib}=[];
   for it = 1:dim(3)
      vtmp = squeeze(var(:,:,it));
      vtmp(iinocn{ib}) = NaN; % set non-basin indices' values to NaN
      iivar{ib}   = cat(1,iivar{ib},find(vtmp>0)+(it-1)*ngp); % search the nonzero indices in the basin for deep convective precip
   end
end   
%save('preccdzm_index.mat','iipzm','itimes'); % each iipzm{ib} can index the 3D variable as var{idv.pzm}(iipzm{ib}), no need to specify itimes
