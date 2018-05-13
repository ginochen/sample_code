% var = loadvar(varname, Case, factor, minusNtime,start,count,stride)
%
% 850mb, ilev=24
% 200mb, ilev=14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D var 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cam   .rh0. var starts from 0001-01-14-03600 to 0001-01-29-64800
% cam   .r.   var starts from 0001-01-14-03600 to 0001-01-29-64800
% spcam .rh0. var starts from 0001-01-14-03600 to 0001-01-29-63000
% spcam .r.   var starts from 0001-01-14-01800 to 0001-01-29-63000
idv.dt=1;  varname{idv.dt} = 'DTCOND';   var{idv.dt} = loadvar(varname{idv.dt},camCase,spd,1);% (K/day) <---- only moist parametrized heating
idv.spdt=30; varname{idv.spdt} = 'SPDT';     var{idv.spdt} = loadvar(varname{idv.spdt},spcamCase,spd,0);% T tendency due to CRM (without QRL+QRS radiative heating)
%for it=1:numel(itimes)
%   dtcond_pm(:,:,it)  = pwgtave(var{idv.dt}(:,:,:,itimes(it)),iilev,dlev);
%   spdt_pm(:,:,it)  = pwgtave(var{idv.spdt}(:,:,:,itimes(it)),iilev,dlev);
%end
%iv=15; varname{iv} = 'OMEGAT';   var{iv} = loadvar(varname{iv},camCase,1,1);%  (3D vertical heat flux)
% ZM mass flux at the 850 to 985 pressure 'interface' a total of 31 ilevs
idv.cmfmcdzm=6; varname{idv.cmfmcdzm} = 'CMFMCDZM'; var{idv.cmfmcdzm} = loadvar(varname{idv.cmfmcdzm},camCase,1,1); % just read level from 850 to 985 skip 1000
%idv.cmfmcdzm=6; varname{idv.cmfmcdzm} = 'CMFMCDZM'; cmfmcdzm_pi24pi30 = loadvar(varname{idv.cmfmcdzm},camCase,1,1,[1,1,24,1],[Inf,Inf,7,Inf]); % just read level from 850 to 985 skip 1000
idv.cmfmcdzm=6; varname{idv.cmfmcdzm} = 'CMFMCDZM'; cmfmcdzm_pi14pi24 = loadvar(varname{idv.cmfmcdzm},camCase,1,1,[1,1,14,1],[Inf,Inf,11,Inf]); % just read level from 850 to 985 skip 1000
idv.cmfmc=7; varname{idv.cmfmc} = 'CMFMC'; var{idv.cmfmc} = loadvar(varname{idv.cmfmc},camCase,1,1); % just read level from 850 to 985 skip 1000
idv.cmfmc=7; varname{idv.cmfmc} = 'CMFMC'; cmfmc_pi24pi30 = loadvar(varname{idv.cmfmc},camCase,1,1,[1,1,24,1],[Inf,Inf,7,Inf]); % just read level from 850 to 985 skip 1000
%idv.cmfmc=7; varname{idv.cmfmc} = 'CMFMC'; cmfmc_pi14pi24 = loadvar(varname{idv.cmfmc},camCase,1,1,[1,1,14,1],[Inf,Inf,11,Inf]); % just read level from 850 to 985 skip 1000
% vertically summed mass flux
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/p-ave/cmfmcdzm_cmfmc_relhum_p-ave.mat')
%cmfmcdzm_pi24pi30m(:,:,:) =  squeeze(mean(cmfmcdzm_pi24pi30,3));%squeeze(mean(var{idv.cmfmcdzm},3));
%cmfmcdzm_pi14pi24m(:,:,:) =  squeeze(mean(cmfmcdzm_pi14pi24,3));%squeeze(mean(var{idv.cmfmcdzm},3));
%cmfmc_pi24pi30m(:,:,:)    =  squeeze(mean(cmfmc_pi24pi30,3));%squeeze(mean(var{idv.cmfmcdzm},3));
%cmfmc_pi21pi30m(:,:,:)    =  squeeze(mean(var{idv.cmfmc}(:,:,21:30,:),3));%squeeze(mean(var{idv.cmfmcdzm},3));
%cmfmc_pi14pi21m(:,:,:)    =  squeeze(mean(var{idv.cmfmc}(:,:,14:21,:),3));%squeeze(mean(var{idv.cmfmcdzm},3));
cmfmcdzm_pi21pi30m(:,:,:)    =  squeeze(mean(var{idv.cmfmcdzm}(:,:,21:30,:),3));%squeeze(mean(var{idv.cmfmcdzm},3));
cmfmcdzm_pi14pi21m(:,:,:)    =  squeeze(mean(var{idv.cmfmcdzm}(:,:,14:21,:),3));%squeeze(mean(var{idv.cmfmcdzm},3));
%cmfmc_pi14pi24m(:,:,:)    =  squeeze(mean(cmfmc_pi14pi24,3));
% save('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/p-ave/cmfmcdzm_cmfmc_relhum_p-ave.mat','relhum_p14p24m','relhum_p24p30m','cmfmc_pi14pi24m','cmfmc_pi24pi30m','cmfmcdzm_pi24pi30m','cmfmcdzm_pi14pi24m')
%
%idv.rh=3; varname{idv.rh} = 'RELHUM';   var{idv.rh} = loadvar(varname{idv.rh},camCase,0.01,1);
idv.rh=3; varname{idv.rh} = 'RELHUM';    relhum_p24p30 = loadvar(varname{idv.rh},camCase,0.01,1,[1,1,24,1],[Inf,Inf,Inf,Inf]);
relhum_p24p30m=squeeze(mean(relhum_p24p30,3));
idv.rh=3; varname{idv.rh} = 'RELHUM';    relhum_p14p24 = loadvar(varname{idv.rh},camCase,0.01,1,[1,1,14,1],[Inf,Inf,10,Inf]);
for it=1:numel(itimes)
   relhum_p24p30m(:,:,it)  = pwgtave(relhum_p24p30(:,:,:,itimes(it)),1:size(relhum_p24p30,3),dlev(24:30));
   relhum_p14p24m(:,:,it)  = pwgtave(relhum_p14p24(:,:,:,itimes(it)),1:size(relhum_p14p24,3),dlev(14:24));
end
%idv.qrs=6;  varname{idv.qrs} = 'QRS';      var{idv.qrs} = loadvar(varname{idv.qrs},spcamCase,spd); % (mm/day)
%idv.qrl=7;  varname{idv.qrl} = 'QRL';      var{idv.qrl} = loadvar(varname{idv.qrl},spcamCase,spd); %(mm/day) <--- is this SPQRL?
%var{2}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.DTCOND.nc'],  'DTCOND')*spd; % (K/day) <---- this is SPDT + SPQRL + SPQRS (QRS_CAM + QRL_CAM ?)
%var{3}  =-ncread([archive      'cam_diff_m2005/m2005.diff.cam.rh0.DTCOND.nc'],'DTCOND')*spd; % spcam-cam (K/day) <--- incorrect
%var{14} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.OMEGA.nc'],              'OMEGA'); %  (3D vertical velocity)
%var{21} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.U.nc'],                  'U');
%var{22} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.V.nc'],                  'V');
%contourf(var{1}(:,:,end,1)-var{6}(:,:,end,1)-var{7}(:,:,end,1)-var{30}(:,:,end,1))
idv.e=31; varname{idv.e} = 'E';        var{idv.e} = loadvar(varname{idv.e},camCase,1,0);% nt=754, b/c var{30} - var{1}(:,:,:,1:end-1); 
idv.ea=32; varname{idv.ea} = 'EABS';     var{idv.ea} = abs(var{idv.e});
for it=701:numel(itimes)
%   e_pm(:,:,it)  = pwgtave(var{31}(:,:,:,itimes(it)),iilev,dlev);
%   eabs_pm(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),iilev,dlev);
%   e_p24p30m(:,:,it)  = pwgtave(var{31}(:,:,:,itimes(it)),[24:30],dlev);
%   eabs_p24p30m(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),[24:30],dlev);
   eabs_p14p21m(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),[14:21],dlev);
   eabs_p21p30m(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),[21:30],dlev);
%   eabs_p24p30m(:,:,it) = pwgtave(var{32}(:,:,:,itimes(it)),[14:24],dlev);
end
%iv=33; varname{iv} = 'AWNC';     var{iv} = loadvar(varname{iv},camCase,1,1);
%iv=34; varname{iv} = 'AWNI';     var{iv} = loadvar(varname{iv},camCase,1,1);
idv.ci=35; varname{idv.ci} = 'CLDICE';   var{idv.ci} = loadvar(varname{idv.ci},camCase,1,1);
idv.clq=36; varname{idv.clq} = 'CLDLIQ';   var{idv.clq} = loadvar(varname{idv.clq},camCase,1,1);
idv.cld=37; varname{idv.cld} = 'CLOUD';    var{idv.cld} = loadvar(varname{idv.cld},camCase,1,1);
idv.ccld=38; varname{idv.ccld} = 'CONCLD';   var{idv.ccld} = loadvar(varname{idv.ccld},camCase,1,1);
%iv=39; varname{iv} = 'FICE';     var{iv} = loadvar(varname{iv},camCase,1,1); % 1's and 0's not interesting
%iv=40; varname{iv} = 'ICIMR';    var{iv} = loadvar(varname{iv},camCase,1,1);
%iv=41; varname{iv} = 'ICLDIWP';  var{iv} = loadvar(varname{iv},camCase,1,1);
idv.icwmr=42; varname{idv.icwmr} = 'ICWMR';    var{idv.icwmr} = loadvar(varname{idv.icwmr},camCase,1,1); % in-cloud water mixing ratio
idv.lccld=43; varname{idv.lccld} = 'LCLOUD';   var{idv.lccld} = loadvar(varname{idv.lccld},camCase,1,1);
idv.ni=44; varname{idv.ni} = 'NUMICE';   var{idv.ni} = loadvar(varname{idv.ni},camCase,1,1);
idv.nl=45; varname{idv.nl} = 'NUMLIQ';   var{idv.nl} = loadvar(varname{idv.nl},camCase,1,1);
idv.spp=46; varname{idv.spp} = 'SPPFLX';   var{idv.spp} = loadvar(varname{idv.spp},spcamCase,spd*1000,0); % Precip flux for CRM
idv.rhcrm=47; varname{idv.rhcrm} = 'RELHUM';   var{idv.rhcrm} = loadvar(varname{idv.rhcrm},spcamCase,0.01,0); % RELHUM for CRM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2D var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the ocean data used in the simulation
% /projects/rsmas/kirtman/cesm_inputdata/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_clim_c040926.nc
% fid = fopen(['/bkirtman4/gchen/cesm_spcam/archive/F_2000_4SPCAM_m200501/ocn/rest' 'F_2000_4SPCAM_m200501.docn.rs1.0001-01-19-00000.bin'],'r');
% SST = fread(fid)
idv.precc=4;     varname{idv.precc} = 'PRECC';       var{idv.precc}   = loadvar(varname{idv.precc},camCase,spd*1000,1);% (mm/day) 
idv.prect=5;     varname{idv.prect} = 'PRECT';       var{idv.prect}   = loadvar(varname{idv.prect},camCase,spd*1000,1);% (mm/day) 
idv.cape=8;      varname{idv.cape} = 'CAPE';         var{idv.cape}    = loadvar(varname{idv.cape},camCase,1,1); 
idv.cin=9;       varname{idv.cin} = 'CIN';           var{idv.cin}     = loadvar(varname{idv.cin},camCase,1,1); 
idv.pzm=10;      varname{idv.pzm} = 'PRECCDZM';      var{idv.pzm}     = loadvar(varname{idv.pzm},camCase,spd*1000,1);% (mm/day) 
idv.psh=11;      varname{idv.psh} = 'PRECSH';        var{idv.psh}     = loadvar(varname{idv.psh},camCase,spd*1000,1);% (mm/day) 
idv.cldh=12;     varname{idv.cldh} = 'CLDHGH';       var{idv.cldh}    = loadvar(varname{idv.cldh},camCase,1,1);
idv.cldl=15;     varname{idv.cldl} = 'CLDLOW';       var{idv.cldl}    = loadvar(varname{idv.cldl},camCase,1,1);
idv.cldm=16;     varname{idv.cldm} = 'CLDMED';       var{idv.cldm}    = loadvar(varname{idv.cldm},camCase,1,1);
idv.cldt=17;     varname{idv.cldt} = 'CLDTOT';       var{idv.cldt}    = loadvar(varname{idv.cldt},camCase,1,1);
idv.solin=18;    varname{idv.solin} = 'SOLIN';       var{idv.solin}   = loadvar(varname{idv.solin},spcamCase,1,0);
idv.tgcldcwp=19; varname{idv.tgcldcwp} = 'TGCLDCWP'; var{idv.tgcldcwp}= loadvar(varname{idv.tgcldcwp},camCase,1,1);
idv.tgcldiwp=20; varname{idv.tgcldiwp} = 'TGCLDIWP'; var{idv.tgcldiwp}= loadvar(varname{idv.tgcldiwp},camCase,1,1);
idv.tgcldlwp=21; varname{idv.tgcldlwp} = 'TGCLDLWP'; var{idv.tgcldlwp}= loadvar(varname{idv.tgcldlwp},camCase,1,1);

nv = numel(find((cellfun('isempty', var))==0)); % number of non-empty vars

%var{5}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.PRECIP.nc'],  'PRECC')*spd*1000; % (mm/day)
%var{6}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.QRS.nc'],     'QRS')*spd;
%var{7}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.QRL.nc'],     'QRL')*spd;
%var{6}  = ncread([camCaseDir   'F_2000_4SPCAM_m200501-cape-cin.nc'],                   'CAPE');
%var{7}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.cape_cin.nc', 'CAPE');
%var{8}  = ncread([camCaseDir   'F_2000_4SPCAM_m200501-cape-cin.nc'],                   'CIN'); 
%var{9}  = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.cape_cin.nc'],'CIN'); 
%var{9}  = ncread([spcamCaseDir
%'spcam_actual_m2005_f09f09_branch.cam.rh0.PRECIP.nc'],                             'PRECCDZM')*spd*1000; % same as CAM PRECCDZM
%var{12} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.PRECIP.nc'],             'PRECL')*spd*1000;
%var{13} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.PRECIP.nc'],             'PRECT')*spd*1000;
%var{16} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.TMQ.nc'],                'TMQ'); %  (Total (vertically integrated) precipitable water)
%var{17} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.ATMEINT.nc'],            'ATMEINT');% (Vertically integrated total atmospheric energy)
%var{19} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.r.U.nc'],                    'U');
%var{20} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.r.V.nc'],                    'V');
%var{23} = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.r.U.nc'],         'U');
%var{24} = ncread([spcamCaseDir 'spcam_actual_m2005_f09f09_branch.cam.rh0.U.nc'],       'U');
%var{25} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.US.nc'],                 'US');
%var{26} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.FSDS.nc'],               'FSDS');
%var{27} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.SOLIN.nc'],              'SOLIN'); % missing data every two days
% solin(i) = tot_irrad*1.e3_r8*eccf*coszrs(i) where i is the day light column index
%var{28} = ncread([camCaseDir   'F_2000_4SPCAM_m200501.cam.rh0.SRFRAD.nc'],             'SRFRAD');
%var{29} = ncread([camCaseDir  'F_2000_4SPCAM_m200501.cam.rh0.CMFDT.nc'],              'CMFDT'); % all zero value, don't use
