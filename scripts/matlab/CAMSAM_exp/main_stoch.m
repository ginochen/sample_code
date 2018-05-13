run ~/scripts/matlab/startup.m
%
archive   = ['/projects/rsmas/kirtman/gchen/archive/'];
%archive   = ['/bkirtman4/gchen/cesm_spcam/archive/'];
camCase.dir   = [archive 'F_2000_4SPCAM_m200501/atm/rest/'];
spcamCase.dir = [archive 'spcam_cam_rh0_m2005/atm/rest/'];
camCase.name = 'F_2000_4SPCAM_m200501.cam.rh0';
spcamCase.name = 'spcam_actual_m2005_f09f09_branch.cam.rh0';
%
%
%run loadparm.m
%run loadvars.m
%run corr_space.m
%run corr_time_map.m
