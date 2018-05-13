function iivar = var_index(var)
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/index/basin_index.mat'); 
% tropics non-zero PRECCDZM index 20S < lat < 20N, basin-wise
ngp = size(var,1)*size(var,2); %number of grid points
for ib = numel(basin)
   iivar{ib}=[];
   for it = 1:numel(itimes)
      vtmp = squeeze(var(:,:,it));
      vtmp(iinocn{ib}) = NaN; % set non-basin indices' values to NaN
      iivar{ib}   = cat(1,iivar{ib},find(vtmp>0)+(it-1)*ngp); % search the nonzero indices in the basin for deep convective precip
   end
end   
%save('preccdzm_index.mat','iipzm','itimes'); % each iipzm{ib} can index the 3D variable as var{idv.pzm}(iipzm{ib}), no need to specify itimes
