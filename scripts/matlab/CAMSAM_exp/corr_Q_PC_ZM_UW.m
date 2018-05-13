% Purpose: Correlation between the PC samples of CAM and SPCAM.
% Caution: only the first three modes are considered
%
% - Q_PC conditioned on UW_shallow only and ZM_deep only
% - concat all basins into one matrix to calc the overall correlation
load '~/archive/matlab/var/var/Q_ZMexcUW_UWexcZM_basinwise-svd.mat'
%load ~/archive/matlab/var/var/Q_ZMincUW_UWexcZM_basinwise-svd_tminus1.mat
%load '~/archive/matlab/var/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
npc_s=[]; npc_sps=[]; 
npc_z=[]; npc_spz=[]; 
pc_sall=[]; pc_zall=[]; 
pc_spsall=[]; pc_spzall=[]; 
for ib=1:4
  for imode=1:3
     [npc_s{ib}(:,imode)   bsign_s(ib,imode)]   = pc_sign_correct(pc_s{ib}(:,imode),  svec_s{ib}(:,imode));
     [npc_sps{ib}(:,imode) bsign_sps(ib,imode)] = pc_sign_correct(pc_sps{ib}(:,imode),svec_sps{ib}(:,imode));
     [npc_z{ib}(:,imode)   bsign_z(ib,imode)]   = pc_sign_correct(pc_z{ib}(:,imode),  svec_z{ib}(:,imode));
     [npc_spz{ib}(:,imode) bsign_spz(ib,imode)] = pc_sign_correct(pc_spz{ib}(:,imode),svec_spz{ib}(:,imode));
  end
  pc_sall   = cat(1,npc_s{ib},  pc_sall); % cat all basins into one matrix 
  pc_spsall = cat(1,npc_sps{ib},pc_spsall); 
  pc_zall   = cat(1,npc_z{ib},  pc_zall); 
  pc_spzall = cat(1,npc_spz{ib},pc_spzall); 
end
[rho_s pval_s]=corr(pc_sall,pc_spsall)
[rho_z pval_z]=corr(pc_zall,pc_spzall)
%rho_s =
%   -0.0493    0.1672   -0.0796
%   -0.0212   -0.0513   -0.0019
%    0.1290   -0.1375    0.0272
%rho_z =
%    0.5400   -0.3208   -0.0673
%   -0.1231    0.0888    0.0061
%   -0.1859    0.2366   -0.0711

function [ pc_new, bsign ] = pc_sign_correct(pc,svec)
   % Correct the sign issue due to svd (flipping the singular vector randomly). 
   % The all-basin correlation between spcam and cam will be very low if the
   % signs in each basin aren't corrected.
   % pc_new : sign-corrected pc
   bsign    = sign(corr(mean(pc)*svec,svec)); 
   pc_new   = pc*bsign;
end   
