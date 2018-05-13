function [Rround P] = mcs_pc_corr(season,vname,dlat,ilatzone,pc,svec,dosave)
%[Rr.dt P.dt] = mcs_pc_corr(season,'dt',dlat,ilatzone,pc,svec,dosave)
%[Rr.dq P.dq] = mcs_pc_corr(season,'dq',dlat,ilatzone,pc,svec,dosave)
ii=1;
latzone = {'trop','midlat','pole','global'};
diro = ['/Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/figure/' dlat '/' season '/' latzone{ilatzone+1} '/']
nmode=2; % EOF3 variance is approx 10% or less for both camdt/camdq or spdt/spdq
eval(sprintf('v1 = pc.sp%s;',vname));
eval(sprintf('v2 = pc.cam%s;',vname));
for ilnd = [0 1]
  for ilife = [1 10 20]
    [R{ii} P{ii}] = corr(v1{ii}(:,1:nmode),v2{ii}(:,1:nmode));
    ii=ii+1;
end; end
Rround = cellfun(@(x) round(x,2),R,'Un',0); % corr rounded to 2nd decimal
if dosave
  eval(sprintf('save([diro ''corr_sp%s_cam%s_pc.mat''],''R'',''Rround'',''P'',''pc'',''svec'');',vname,vname))
  disp('saved')
end
