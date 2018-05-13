
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-correlated, time-standard deviated 3D map    
% (zero indices excluded according to vtmp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
do_pwgave=0;
iilev=[1:29];
ivv=35;
%ivv = [42,43,44,45];
itimes=[1:700];
fileE='/projects/rsmas/kirtman/gchen/archive/matlab/figure/rho_map/E/rho_p_map_e_30lev_test.mat';
fileEABS='/projects/rsmas/kirtman/gchen/archive/matlab/figure/rho_map/EABS/rho_p_map_eabs_30lev_test.mat';
if exist(fileE,'file') | exist(fileEABS,'file')
   error('fileE exist')
end
save(fileE,'ilat','ilon','iilev','itimes')
save(fileEABS,'ilat','ilon','iilev','itimes')
for iv=1:numel(ivv)
   tic
   if ( do_pwgave )
      for it=1:numel(itimes)
         vtmp(:,:,it)=pwgtave(var{ivv(iv)}(:,:,:,itimes(it)),iilev,dlev);
      end
   else % ndim(var{iv})==3
      if ( ndims(var{ivv(iv)})==3 )
         vtmp = var{ivv(iv)}(:,:,itimes);
      else %ndims == 4
         vtmp = var{ivv(iv)}(:,:,:,itimes);
      end
   end
   [rhomap_e{iv} pmap_e{iv}] = corrcoef_time_map_nzero(var{31}(:,:,:,itimes),vtmp,ilat,ilon,iilev);      
   [rhomap_eabs{iv} pmap_eabs{iv}] = corrcoef_time_map_nzero(var{32}(:,:,:,itimes),vtmp,ilat,ilon,iilev);      
   toc
   vname{iv} = varname{ivv(iv)};
   save(fileE,'rhomap_e','pmap_e','vname','-append')
   save(fileEABS,'rhomap_eabs','pmap_eabs','vname','-append')
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-correlated, time-standard deviated 2D map or 3 level map   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ivv = [35,36,38,42]
illev{1}=[1:18]; illev{2}=[18:24]; illev{3}=[24:30];% [1,18,24,30] % [2 400 850 1000] mb
fileE='/projects/rsmas/kirtman/gchen/archive/matlab/figure/rho_map/E/rho_p_map_cld_e_3lev.mat',
fileEABS='/projects/rsmas/kirtman/gchen/archive/matlab/figure/rho_map/EABS/rho_p_map_cld_eabs_3lev.mat',
save(fileE,'ilat','ilon','illev','itimes')
save(fileEABS,'ilat','ilon','illev','itimes')
% convert E and EABS into p-weighted 2D maps
for ill = 1:numel(illev) % split into three levels and do correlation
   for it=1:numel(itimes)
      etmp{ill}(:,:,it)  = pwgtave(var{31}(:,:,:,itimes(it)),     ilev(illev{ill}),dlev(illev{ill}));
      eatmp{ill}(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),     ilev(illev{ill}),dlev(illev{ill}));
   end
end
% threee level maps
for iv=1:numel(ivv)
   tic
   for ill = 1:numel(illev) % split into three levels and do correlation
      for it=1:numel(itimes)
         vtmp(:,:,it)  = pwgtave(var{ivv(iv)}(:,:,:,itimes(it)),ilev(illev{ill}),dlev(illev{ill}));
      end
      [rhomap_e{iv}(:,:,ill)    pmap_e{iv}(:,:,ill)]    = corrcoef_time_map_nzero(etmp{ill}, vtmp,ilat,ilon);      
      [rhomap_eabs{iv}(:,:,ill) pmap_eabs{iv}(:,:,ill)] = corrcoef_time_map_nzero(eatmp{ill},vtmp,ilat,ilon);      
   end
   vname{iv} = varname{ivv(iv)};
   save(fileE,'rhomap_e','pmap_e','vname','-append')
   save(fileEABS,'rhomap_eabs','pmap_eabs','vname','-append')
   toc
end
