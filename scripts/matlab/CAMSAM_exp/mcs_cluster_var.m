function mcs_cluster_var(fci,fin,fout,fmcs,fparm,time,season)
%function mcs_cluster_var(fin,season)
% Purpose: obtain the variables for the extra lag-lead index for clusters, 
%   save each time in separate files. Use mcs_cluster_combineExt.m to combine
%   these files later.
donew = 0;
  disp('make sure to uncomment after doing new variables')
  close all force
  run ~/scripts/matlab/startup.m
%  time        = fin(end-18:end-3) % get the time of nc file
%  spcase      = 'F_2000_SPCAM_m2005_3hrly1';
%  spcase      = 'F_2000_SPCAM_m2005_3hrly_f09f09_branch';
%  spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
%  spArchive   = ['/glade2/scratch2/ginochen/archive/matlab/'];
%  spHist.dir  = [spArchive spcase '/atm/hist/']
  %season = {'MAM','JJA','SON','DJF'};
%  fout = [spHist.dir season '/mcs_cluster_var/mcs_cluster.' time '.mat']; disp(fout)
  load(fmcs,'t','t4Cl','nt4Cl','mcsn','mcsnll','nCl4nt','nCl','nt4ClRowInd','mcsillt4Cl','mcsllx4Cl');
  it = find(strcmp(time,t))
  if ~mcsn(it); error('no mcs in this time, exiting'); end
  nt = numel(t);
  if ismember(it,[1,nt]); error(['time index = ' num2str(it) '. make sure it starts no less than 2 and end before ' num2str(nt)]); end
%  fin = [spHist.dir spcase '.cam.h1.' t{1} '.nc']; 
  parm = mcs_cluster_parm(fin);
  if ~exist(fparm,'file')
    save(fparm,'parm');
  end
  if ~exist(fout,'file') 
    save(fout,'parm'); 
    vo = [];
  else
    load(fout,'vo','parm')
  end
%{
  vo.B = [];
%}  
tic
%  for it=2:numel(t)-1
    fin0 = [fci '.' t{it} '.nc'];
    fin1 = [fci '.' t{it-1} '.nc'] ;
    fin2 = [fci '.' t{it+1} '.nc'] ; 
    [vi0] = readallvar(fin0,donew);
    [vi1] = readallvar(fin1,donew);
    [vi2] = readallvar(fin2,donew);
    for ic = 1:nCl
%      if any(ismember(t4Cl(ic,:),[1 nt])) | nt4Cl(ic)==1; continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
      if any(ismember(t4Cl(ic,:),[1 nt])); iCends(ic)=1; continue; else iCends(ic)=0; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
      if (it == mcsillt4Cl{ic}{1}(1,3)) % cluster's starting time index from the (1,3)'s 1 represents the 1st subcluster
        for ill=1:mcsnll(ic,1) % go through the ilon-ilat at iit time for cluster ic
          [vo] = calcallvar(1,ic,ill,parm,mcsillt4Cl{ic}{1},mcsllx4Cl{ic}{1}, vi1, vo, donew); 
          [vo] = calcallvar(2,ic,ill,parm,mcsillt4Cl{ic}{1},mcsllx4Cl{ic}{1}, vi0, vo, donew);
          if (it == mcsillt4Cl{ic}{end}(1,3)) % for short lived clusters with only one time index 
            [vo] = calcallvar(3,ic,ill,parm,mcsillt4Cl{ic}{end},mcsllx4Cl{ic}{end}, vi2, vo, donew);
          end
        end
      elseif (it < mcsillt4Cl{ic}{end}(1,3) & it > mcsillt4Cl{ic}{1}(1,3))% it's in between start and end of a cluster
        itt = find(t4Cl(ic,:)==it); % t4Cl contains time indices evolution of ic-th cluster from the set 1:nt
        for ill=1:mcsnll(ic,itt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
          [vo] = calcallvar(itt+1,ic,ill,parm,mcsillt4Cl{ic}{itt},mcsllx4Cl{ic}{itt}, vi0, vo, donew);
        end
      elseif (it == mcsillt4Cl{ic}{end}(1,3)) % cluster's ending time index
        ntt = nt4Cl(ic);
        for ill=1:mcsnll(ic,ntt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
          [vo] = calcallvar(ntt+1,ic,ill,parm,mcsillt4Cl{ic}{end},mcsllx4Cl{ic}{end}, vi0, vo, donew); 
          [vo] = calcallvar(ntt+2,ic,ill,parm,mcsillt4Cl{ic}{end},mcsllx4Cl{ic}{end}, vi2, vo, donew);
        end
      end % if loop of calc variables 
    end % for loop of cluster indices
    save(fout, 'vo','-append','-v6')
    save(fmcs, 'iCends','-append')
    close all force
toc
%  end % for loop of time indices



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [vi] = readallvar(fin, donew)
if donew 
%
else
  vi.u_cam    = flipdim(ncread(fin,'U'),3);
  vi.v_cam    = flipdim(ncread(fin,'V'),3);
  vi.u_crm    = squeeze(ncread(fin,'CRM_U')); % (lon, lat, nx, ny, nz, nt)
  vi.phis     = squeeze(ncread(fin,'PHIS')); 
  vi.spdt     = flipdim(ncread(fin,'SPDT'),3); % sp mass flux [kg/m^2/s] DTCOND (includes QRS + QRL and some weird offsets, not just moist process)
  vi.spdq     = flipdim(ncread(fin,'SPDQ'),3); % sp mass flux [kg/m^2/s] no difference with DCQ
%  vi.spmc     = flipdim(ncread(fin,'SPMC'),3); % sp mass flux [kg/m^2/s]
%  vi.spmcup   = flipdim(ncread(fin,'SPMCUP'),3); % sp updraft mass flux 
%  vi.spmcdn   = flipdim(ncread(fin,'SPMCDN'),3); % sp downdraft mass flux 
%  vi.spmc     = flipdim(ncread(fin,'SPMC'),3); % sp mass flux [kg/m^2/s]
  vi.Ts       = squeeze(ncread(fin,'TS'));
% for doing LI
%  vi.dtcond = flipdim(ncread(fin,'DTCOND'),3); % "T tendency - moist processes" [K/s]
%  vi.dcq    = flipdim(ncread(fin,'DCQ'),3); % "Q tendency due to moist processes" [kg/kg/s]
  vi.flnt   = squeeze(ncread(fin,'FLNT')); % "Net longwave flux at top of model" [W/m^2]
  vi.fsnt   = squeeze(ncread(fin,'FSNT')); % "Net solar flux at top of model" [W/m^2]

  vi.T_crm    = squeeze(ncread(fin,'CRM_T')); % (lon, lat, nx, ny, nz, nt)
  vi.prec_crm = squeeze(ncread(fin,'CRM_PREC'))*3.6e6; %[mm/hr]
  vi.qv_crm       = squeeze(ncread(fin,'CRM_QV')); % (lon, lat, nx, ny, nz, nt)
  vi.qc_crm       = squeeze(ncread(fin,'CRM_QC')); % (lon, lat, nx, ny, nz, nt)
  vi.qi_crm       = squeeze(ncread(fin,'CRM_QI')); % (lon, lat, nx, ny, nz, nt)
  vi.qr_crm       = squeeze(ncread(fin,'CRM_QR')); % (lon, lat, nx, ny, nz, nt)
  vi.w_crm    = squeeze(ncread(fin,'CRM_W')); % surface w_crm = 0, all zeros
  %vi.camdt      = vi.cmfdt+vi.zmdt+vi.zmmtt+(vi.mpdt+vi.macpdt)/1004; % cam parameterized heating, the original DTCOND
  vi.camdt    = flipdim(ncread(fin,'ZMDT')+ncread(fin,'EVAPTZM')+ncread(fin,'CMFDT')+(ncread(fin,'MPDT')+ncread(fin,'MACPDT'))/1004, 3);
  vi.camdq    = flipdim(ncread(fin,'ZMDQ')+ncread(fin,'EVAPQZM')+ncread(fin,'CMFDQ')+(ncread(fin,'MPDQ')+ncread(fin,'MACPDQ')), 3);
%  vi.precc_cam  = ncread(fin,'PRECC')*3.6e6; % convective prec rate [m/s]->[mm/hr] <<<<<<<<< PRECC = mean(CRM_PREC) no need to save here
%  vi.precl_cam  = ncread(fin,'PRECL')*3.6e6; % stratiform prec rate [m/s]->[mm/hr]
%  vi.precsc  = ncread(fin,'PRECSH')*3.6e6; % stratiform prec rate is PRECSED, but not saved. [m/s]->[mm/hr]
%  vi.cmfmcdzm = flipdim(ncread(fin,'CMFMCDZM'),3); % [kg/m^2/s] 31 levels at the interfaces
%  vi.cmfmc    = flipdim(ncread(fin,'CMFMC'),3); % [kg/m^2/s] 31 levels at the interfaces
  vi.zmdt     = flipdim(ncread(fin,'ZMDT'),3); % cam parameterized heating, the original DTCOND
  vi.evaptzm  = flipdim(ncread(fin,'EVAPTZM'),3);
  vi.cmfdt    = flipdim(ncread(fin,'CMFDT'),3);
  vi.mpdt     = flipdim(ncread(fin,'MPDT'),3)/1004;
  vi.macpdt   = flipdim(ncread(fin,'MACPDT'),3)/1004;
  vi.zmdq     = flipdim(ncread(fin,'ZMDQ'),3); % cam parameterized heating, the original DTCOND
  vi.evapqzm  = flipdim(ncread(fin,'EVAPQZM'),3);
  vi.cmfdq    = flipdim(ncread(fin,'CMFDQ'),3);
  vi.mpdq     = flipdim(ncread(fin,'MPDQ'),3);
  vi.macpdq   = flipdim(ncread(fin,'MACPDQ'),3);
%  vi.precz    = ncread(fin,'PRECZ')*3.6e6; %  [m/s]->[mm/hr]

%  vi.v        = flipdim(ncread(fin,'V'),3);
  vi.ps       = squeeze(ncread(fin,'PS')); 
  vi.z3       = flipdim(ncread(fin,'Z3'),3);
  vi.T        = flipdim(ncread(fin,'T'),3);
  vi.q        = flipdim(squeeze(ncread(fin,'Q')),3); 
  vi.qt        = flipdim(squeeze(ncread(fin,'QT')),3); 
%  vi.cldliq   = flipdim(squeeze(ncread(fin,'CLDLIQ')),3); 
%  vi.cldice   = flipdim(squeeze(ncread(fin,'CLDICE')),3); 
%  vi.qrain    = flipdim(squeeze(ncread(fin,'QRAIN')),3); 
  vi.omega    = flipdim(squeeze(ncread(fin,'OMEGA')),3); 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vo] = calcallvar(iitt,ic,ill,parm,mcsillt4Cl,mcsllx4Cl, vi, vo, donew) 
%  z0 = vi.phis(l1,l2)/parm.g; % surface geopotential height
%  z30 = squeeze(vi.z3(l1,l2,1:parm.nz))-z0; % geopotential height w.r.t local surface height
  l1 = mcsillt4Cl(ill,1); % lon
  l10 = l1-1; % used for calculating 1st deriv
  l11 = l1+1; 
  if l1 == length(parm.lon) % at right bd
    l11 = 1;
  elseif l1 == 1
    l10 = length(parm.lon); % at left bd
  end
  l2 = mcsillt4Cl(ill,2); % lat
  l20 = l2-1;
  l21 = l2+1;
  vo.ilon{ic}{iitt}(ill) = l1;
  vo.ilat{ic}{iitt}(ill) = l2;
if donew
disp('uncomment [ctype, u_crm, u_cam after rerunning')
%  if exist('vo.u_cami','var'); error('variable exists, exiting...'); end
%  vo.u_cam{ic}{iitt}(:,ill) = vi.u_cam(l1,l2,1:parm.nz);
%  vo.du03km_cam{ic}{iitt}(ill) = u_cam(parm.i3km)-vi.u_cam(l1,l2,1);
%  vo.du03km_cam_abs{ic}{iitt}(ill) = abs(vo.du03km_cam{ic}{iitt}(ill));
else
  [pi, pm, dpi, dpm]  = hybrid2p(parm.p0, vi.ps(l1,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
%  qv_crmLS = squeeze(mean(vi.qv,3)); % large-scale qv vapor mixing ratio for crm 
%  T_crmLS  = squeeze(mean(vi.T_crm,3)); % large-scale qv vapor mixing ratio for crm 
  for iz = 1:parm.nz
%    relhum_crm(iz) = specific2relhum(T_crmLS(l1,l2,iz),pm(iz),qv_crmLS(l1,l2,iz)/(1+qv_crmLS(l1,l2,iz)));
%    thetae_crm(iz) = gettheta_e(pm(iz),T_crmLS(l1,l2,iz),qv_crmLS(l1,l2,iz));
    rho1d(iz)  = density_temp(squeeze(vi.T(l1,l2,iz)),pm(iz),vi.q(l1,l2,iz)./(1-vi.q(l1,l2,iz)),[],[],squeeze(vi.qt(l1,l2,iz)), 2); % midpoint density, disregard the rain effect here
    for ix = 1:parm.nx
      [rho2d(ix,iz), T_rho(ix,iz), rh_crm(ix,iz)]  = density_temp(squeeze(vi.T_crm(l1,l2,ix,iz)), pm(iz),...
          squeeze(vi.qv_crm(l1,l2,ix,iz)),squeeze(vi.qc_crm(l1,l2,ix,iz)),squeeze(vi.qr_crm(l1,l2,ix,iz)),[], 1); % midpoint density
      th_rho(ix,iz) = pottemp(T_rho(ix,iz),pm(iz));
%      if iz==1 % zonal averaged vertical div equals zero since zonal averaged div equals zero
%        div_crm(ix,1) = 9.8*(rho2d(ix,iz)*vi.w_crm(l1,l2,ix,iz)-0)/dpm(iz); % surface w=0, the negative sign is already multiplied in dpm
%      else
%        div_crm(ix,iz) = 9.8*(rho2d(ix,iz)*vi.w_crm(l1,l2,ix,iz)-rho2d(ix,iz-1)*vi.w_crm(l1,l2,ix,iz-1))/dpm(iz); 
%      end
    end
  end
  [z1d, dzm1d, dzi1d] = p2z(rho1d,dpi,dpm,parm.g,parm.nz);
  u_cam = interp1(z1d, squeeze(vi.u_cam(l1,l2,1:parm.nz)), parm.zint','linear');
  v_cam = interp1(z1d, squeeze(vi.v_cam(l1,l2,1:parm.nz)), parm.zint','linear');
  vo.u_cami{ic}{iitt}(:,ill) = u_cam;
  vo.v_cami{ic}{iitt}(:,ill) = v_cam;
  th_rho_ave=mean(th_rho,1); % LS ave
  for ix = 1:parm.nx
    [z2d(ix,:), dzm2d(ix,:), dzi2d(ix,:)] = p2z(rho2d(ix,:),dpi,dpm,parm.g,parm.nz);
    [B(ix) buoy(ix,:)] = cold_pool_intensity(z2d(ix,:)',th_rho(ix,:),th_rho_ave,parm);
    LI(ix) = liftedIndex(vi.Ts(l1,l2),vi.ps(l1,l2),vi.qv_crm(l1,l2,ix,1),pm,squeeze(vi.T(l1,l2,:)));
    u_crm(ix,:) = interp1(z2d(ix,:), squeeze(vi.u_crm(l1,l2,ix,1:parm.nz)), parm.zint','linear');
  end
%%%%%%%%%%%%%%% OUTPUT VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [ctype, prate] = conv_strat_separation(squeeze(vi.prec_crm(l1,l2,:)), squeeze(vi.w_crm(l1,l2,:,:)), squeeze(vi.qc_crm(l1,l2,:,:)), squeeze(vi.T_crm(l1,l2,:,:)),pm); 
  vo.frac_precl{ic}{iitt}(ill) = nansum(prate(intersect(find(ctype==0),mcsllx4Cl{ill})))/nansum(prate(mcsllx4Cl{ill}));
  vo.frac_prec{ic}{iitt}(ill,1)= nansum(prate(intersect(find(prate>=0.1),mcsllx4Cl{ill})))/nansum(prate(mcsllx4Cl{ill}));
  vo.frac_prec{ic}{iitt}(ill,2)= nansum(prate(intersect(find(prate>=3),mcsllx4Cl{ill})))/nansum(prate(mcsllx4Cl{ill}));
  vo.frac_prec{ic}{iitt}(ill,3)= nansum(prate(intersect(find(prate>=5),mcsllx4Cl{ill})))/nansum(prate(mcsllx4Cl{ill}));
  vo.frac_precc{ic}{iitt}(ill) = nansum(prate(intersect(find(ctype==1),mcsllx4Cl{ill})))/nansum(prate(mcsllx4Cl{ill}));
  vo.freq_precc{ic}{iitt}(ill) = nansum(ctype==1)/parm.nx;
  vo.freq_precl{ic}{iitt}(ill) = nansum(ctype==0)/parm.nx;
  vo.precc{ic}{iitt}(ill)      = nansum(prate(ctype==1))/sum(ctype==1);
  vo.precl{ic}{iitt}(ill)      = nansum(prate(ctype==0))/sum(ctype==0);

  vo.v_cam{ic}{iitt}(:,ill) = vi.v_cam(l1,l2,1:parm.nz);
  vo.u_cam{ic}{iitt}(:,ill) = vi.u_cam(l1,l2,1:parm.nz);
  vo.u_cami{ic}{iitt}(:,ill) = interp1(z1d, squeeze(vi.u_cam(l1,l2,1:parm.nz)), parm.zint','linear');
  vo.v_cami{ic}{iitt}(:,ill) = interp1(z1d, squeeze(vi.v_cam(l1,l2,1:parm.nz)), parm.zint','linear');
  vo.buoy{ic}{iitt}(:,:,ill)   = buoy; % buoyancy 2D profile
  vo.B{ic}{iitt}(:,ill)        = B;
  vo.Bave{ic}{iitt}(ill)       = mean(B(mcsllx4Cl{ill})); % buoyancy strength, originally vo.B{ic}{iitt}(:,ill) <- not much use saving this
  vo.Bmax{ic}{iitt}(ill)       = max(B(mcsllx4Cl{ill})); % buoyancy strength, originally vo.B{ic}{iitt}(:,ill) <- not much use saving this
  vo.LIave{ic}{iitt}(ill)      = mean(LI(mcsllx4Cl{ill}));
  vo.LImax{ic}{iitt}(ill)      = max(LI(mcsllx4Cl{ill}));
  vo.LImin{ic}{iitt}(ill)      = min(LI(mcsllx4Cl{ill}));
  vo.LI{ic}{iitt}(:,ill)       = LI;
  vo.du03km_cam{ic}{iitt}(ill) = u_cam(parm.i3km)-vi.u_cam(l1,l2,1);
  vo.dv03km_cam{ic}{iitt}(ill) = v_cam(parm.i3km)-vi.v_cam(l1,l2,1);
  vo.dU03km_cam{ic}{iitt}(ill) = sqrt(vo.du03km_cam{ic}{iitt}(ill)^2+vo.dv03km_cam{ic}{iitt}(ill)^2);
  vo.du03km{ic}{iitt}(ill)     = mean(u_crm(:,parm.i3km))-mean(vi.u_crm(l1,l2,:,1));
%  vo.du03km_cam{ic}{iitt}(ill) = u_cam(parm.i3km)-vi.u_cam(l1,l2,1);
  vo.dzm2d{ic}{iitt}(:,:,ill)  = dzm2d;
  vo.dzi2d{ic}{iitt}(:,:,ill)  = dzi2d;
  vo.dzi1d{ic}{iitt}(:,ill)    = dzi1d;
  vo.camdt{ic}{iitt}(:,ill)    = vi.camdt(l1,l2,1:parm.nz);
  vo.camdq{ic}{iitt}(:,ill)    = vi.camdq(l1,l2,1:parm.nz);
  vo.cmfdt{ic}{iitt}(:,ill)    = vi.cmfdt(l1,l2,1:parm.nz);
  vo.cmfdq{ic}{iitt}(:,ill)    = vi.cmfdq(l1,l2,1:parm.nz);
  vo.ctype{ic}{iitt}(:,ill)    = ctype;
  vo.evapqzm{ic}{iitt}(:,ill)  = vi.evapqzm(l1,l2,1:parm.nz);
  vo.evaptzm{ic}{iitt}(:,ill)  = vi.evaptzm(l1,l2,1:parm.nz);
  vo.flnt{ic}{iitt}(ill)       = vi.flnt(l1,l2);
  if isnan(vo.frac_precc{ic}{iitt}(ill)) % no rain in the growth (start -1 time index) or decay (last + 1 time index) phase 
    disp(['there is nan (0/0=NaN) in frac_precc!! maybe just has no rain in the ' num2str(iitt) 'th time in the cluster with total time' ])
  end 
  vo.fsnt{ic}{iitt}(ill)       = vi.fsnt(l1,l2);
%  vo.LI{ic}{iitt}(ill)         = liftedIndex(vi.Ts(l1,l2),vi.ps(l1,l2),max(vi.qv_crm(l1,l2,:,1)),pm,squeeze(vi.T(l1,l2,:)));
  vo.LI{ic}{iitt}(ill)         = liftedIndex(vi.Ts(l1,l2),vi.ps(l1,l2),max(vi.qv_crm(l1,l2,mcsllx4Cl{ill},1)),pm,squeeze(vi.T(l1,l2,:)));
  vo.macpdt{ic}{iitt}(:,ill)   = vi.macpdt(l1,l2,1:parm.nz);
  vo.macpdq{ic}{iitt}(:,ill)   = vi.macpdq(l1,l2,1:parm.nz);
  vo.mcsllx4Cl{ic}{iitt}{ill} = mcsllx4Cl{ill};
  vo.mcsillt4Cl{ic}{iitt}(ill,1:3) = mcsillt4Cl(ill,1:3);
  vo.mpdt{ic}{iitt}(:,ill)     = vi.mpdt(l1,l2,1:parm.nz);
  vo.mpdq{ic}{iitt}(:,ill)     = vi.mpdq(l1,l2,1:parm.nz);
  vo.prate{ic}{iitt}(:,ill)    = prate;
  vo.prect{ic}{iitt}(ill)      = sum(prate)/parm.nx;
  vo.phis{ic}{iitt}(ill)       = vi.phis(l1,l2);
  vo.rho2d{ic}{iitt}(:,:,ill)  = rho2d;
  vo.rho1d{ic}{iitt}(:,ill)    = rho1d;
  vo.spdt{ic}{iitt}(:,ill)     = vi.spdt(l1,l2,1:parm.nz);
  vo.spdq{ic}{iitt}(:,ill)     = vi.spdq(l1,l2,1:parm.nz);
  vo.qv_crm{ic}{iitt}(:,:,ill) = vi.qv_crm(l1,l2,:,1:parm.nz);
  vo.qc_crm{ic}{iitt}(:,:,ill) = vi.qc_crm(l1,l2,:,1:parm.nz);
  vo.qi_crm{ic}{iitt}(:,:,ill) = vi.qi_crm(l1,l2,:,1:parm.nz);
  vo.qr_crm{ic}{iitt}(:,:,ill) = vi.qr_crm(l1,l2,:,1:parm.nz);
  vo.T_crm{ic}{iitt}(:,:,ill)  = vi.T_crm(l1,l2,:,1:parm.nz);
  vo.w_crm{ic}{iitt}(:,:,ill)  = vi.w_crm(l1,l2,:,1:parm.nz);
  vo.z1d{ic}{iitt}(:,ill)      = z1d;
  vo.z2d{ic}{iitt}(:,:,ill)    = z2d;
  vo.z3{ic}{iitt}(:,ill)       = vi.z3(l1,l2,1:parm.nz);
  vo.zmdt{ic}{iitt}(:,ill)     = vi.zmdt(l1,l2,1:parm.nz);
  vo.zmdq{ic}{iitt}(:,ill)     = vi.zmdq(l1,l2,1:parm.nz);
  vo.omega{ic}{iitt}(:,ill)     = vi.omega(l1,l2,1:parm.nz);
end





%  vo.precc_cam{ic}{iitt}(ill)  = vi.precc_cam(l1,l2);
%  vo.precl_cam{ic}{iitt}(ill)  = vi.precl_cam(l1,l2);
%{
  vn2D = {'u','w','T','qv','qc','qi','qr','rh','div','rho','buoy'};
  i2D = containers.Map(vn2D,[1:numel(vn2D)]);
  vn1D = {'u_cam','v_cam','spdt','spdq','camdt','camdq','spmc','spmcup','thetae','div','relhum'};
  i1D = containers.Map(vn1D,[1:numel(vn1D)]);
  vo.rho{ic}{iitt}(:,ill) = interp1(z1d, rho1d', parm.zint','linear'); % forgot to use this for pv
  vo.w{ic}{iitt}(:,ill) = interp1(z1d, -squeeze(vi.omega(l1,l2,1:parm.nz))./(rho1d'*parm.g), parm.zint','linear'); % omega to w
  vo.theta{ic}{iitt}(:,ill) = interp1(z1d,[squeeze(pottemp(squeeze(vi.T(l1,l2,1:parm.nz))',pm))], parm.zint','linear');
  vo.thetaz{ic}{iitt}(:,ill) = (vo.theta{ic}{iitt}(2:end,ill)-vo.theta{ic}{iitt}(1:end-1,ill))/parm.dzi;
  theta10 = interp1(z1d,[squeeze(pottemp(squeeze(vi.T(l10,l2,1:parm.nz))',pm))], parm.zint','linear');
  theta11 = interp1(z1d,[squeeze(pottemp(squeeze(vi.T(l11,l2,1:parm.nz))',pm))], parm.zint','linear');
  vo.thetax{ic}{iitt}(:,ill) = (theta11 - theta10)/parm.dlon(l2);
  theta20 = interp1(z1d,[squeeze(pottemp(squeeze(vi.T(l1,l20,1:parm.nz))',pm))], parm.zint','linear');
  theta21 = interp1(z1d,[squeeze(pottemp(squeeze(vi.T(l1,l21,1:parm.nz))',pm))], parm.zint','linear');
  vo.thetay{ic}{iitt}(:,ill) = (theta21 - theta20)/parm.dlat;
  vo.u_cam{ic}{iitt}(:,ill) = interp1(z1d,[squeeze(vi.u(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.uz{ic}{iitt}(:,ill) = (vo.u_cam{ic}{iitt}(2:end,ill)-vo.u_cam{ic}{iitt}(1:end-1,ill))/parm.dzi;
  u20 = interp1(z1d,[squeeze(vi.u(l1,l20,1:parm.nz))], parm.zint','linear');
  u21 = interp1(z1d,[squeeze(vi.u(l1,l21,1:parm.nz))], parm.zint','linear');
  vo.uy{ic}{iitt}(:,ill) = (u21 - u20)/parm.dlat;
  vo.v_cam{ic}{iitt}(:,ill) = interp1(z1d,[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.vz{ic}{iitt}(:,ill) = (vo.v_cam{ic}{iitt}(2:end,ill)-vo.v_cam{ic}{iitt}(1:end-1,ill))/parm.dzi;
  v10 = interp1(z1d,[squeeze(vi.v(l10,l2,1:parm.nz))], parm.zint','linear');
  v11 = interp1(z1d,[squeeze(vi.v(l11,l2,1:parm.nz))], parm.zint','linear');
  vo.vx{ic}{iitt}(:,ill) = (v11 - v10)/parm.dlon(l2);
  vo.pv{ic}{iitt}(:,ill) = (-vo.vz{ic}{iitt}(:,ill).*vo.thetax{ic}{iitt}(1:end-1,ill)+vo.uz{ic}{iitt}(:,ill).*vo.thetay{ic}{iitt}(1:end-1,ill)+...
          (parm.co(l2) + vo.vx{ic}{iitt}(1:end-1,ill)-vo.uy{ic}{iitt}(1:end-1,ill)).*vo.thetaz{ic}{iitt}(:,ill));

  vo.precz{ic}{iitt}(ill) = vi.precz(l1,l2);
%  vo.cmfmcdzm= interp1(z1d,[squeeze(vi.cmfmcdzm(l1,l2,1:parm.nz))], parm.zint','linear');
%  vo.cmfmc   = interp1(z1d,[squeeze(vi.cmfmc(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.zmdt{ic}{iitt}(:,ill)    = interp1(z1d,[squeeze(vi.zmdt(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.evaptzm{ic}{iitt}(:,ill) = interp1(z1d,[squeeze(vi.evaptzm(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.cmfdt{ic}{iitt}(:,ill)   = interp1(z1d,[squeeze(vi.cmfdt(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.mpdt{ic}{iitt}(:,ill)    = interp1(z1d,[squeeze(vi.mpdt(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.macpdt{ic}{iitt}(:,ill)  = interp1(z1d,[squeeze(vi.macpdt(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.zmdq{ic}{iitt}(:,ill)    = interp1(z1d,[squeeze(vi.zmdq(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.evapqzm{ic}{iitt}(:,ill) = interp1(z1d,[squeeze(vi.evapqzm(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.cmfdq{ic}{iitt}(:,ill)   = interp1(z1d,[squeeze(vi.cmfdq(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.mpdq{ic}{iitt}(:,ill)    = interp1(z1d,[squeeze(vi.mpdq(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.macpdq{ic}{iitt}(:,ill)  = interp1(z1d,[squeeze(vi.macpdq(l1,l2,1:parm.nz))], parm.zint','linear');
  vo.prect{ic}{iitt}(ill)      = mean(vi.prec_crm(l1,l2,:));

  v2D = zeros(nx,parm.nzi,numel(vn2D));
  v1D = zeros(parm.nzi,numel(vn1D));
  for ix = 1:nx
      v2D(ix,1:parm.nzi,:) = interp1(z2d(ix,:)', [ ...
                                     squeeze(vi.u_crm(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.w_crm(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.T_crm(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.qv(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.qc(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.qi(l1,l2,ix,1:parm.nz)),...
                                     squeeze(vi.qr(l1,l2,ix,1:parm.nz)),...
                                     rh_crm(ix,:)',...
                                     div_crm(ix,:)',...
                                     rho2d(ix,:)',...
                                     buoy(ix,:)'], parm.zint,'linear');
  end
  div = mean(div_crm(mcsllx4Cl{ill},:),1);
  vo.u_crm{ic}{iitt}(:,:,ill)   = v2D(:,:,i2D('u'));
  vo.w_crm{ic}{iitt}(:,:,ill)   = v2D(:,:,i2D('w'));
  vo.T_crm{ic}{iitt}(:,:,ill)   = v2D(:,:,i2D('T'));
  vo.qv_crm{ic}{iitt}(:,:,ill)  = v2D(:,:,i2D('qv'));
  vo.qc_crm{ic}{iitt}(:,:,ill)  = v2D(:,:,i2D('qc'));
  vo.qi_crm{ic}{iitt}(:,:,ill)  = v2D(:,:,i2D('qi'));
  vo.qr_crm{ic}{iitt}(:,:,ill)  = v2D(:,:,i2D('qr'));
  vo.div_crm{ic}{iitt}(:,:,ill) = v2D(:,:,i2D('div')); 
  vo.rh_crm{ic}{iitt}(:,:,ill)  = v2D(:,:,i2D('rh'));
  vo.rho_crm{ic}{iitt}(:,:,ill) = v2D(:,:,i2D('rho'));
  v1D = interp1(z1d,[squeeze(vi.u(l1,l2,1:parm.nz)),...
                     squeeze(vi.v(l1,l2,1:parm.nz)),...
                     squeeze(vi.spdt(l1,l2,1:parm.nz)),...
                     squeeze(vi.spdq(l1,l2,1:parm.nz)),...
                     squeeze(vi.camdt(l1,l2,1:parm.nz)),...
                     squeeze(vi.camdq(l1,l2,1:parm.nz)),...
                     squeeze(vi.spmc(l1,l2,1:parm.nz)),...
                     squeeze(vi.spmcup(l1,l2,1:parm.nz)),...
                     thetae_crm',...
                     div',...
                     relhum_crm',...
                     ], parm.zint','linear');
%                                                   div',...        
                                                   %squeeze(vi.spmcdn(l1,l2,1:parm.nz)),...
  %for i=1:numel(vn1D)
  %  eval(sprintf('vo.%s{ic}{iitt}(:,ill) = v1D(:,i1D(''%s''));',vn1D{i},vn1D{i}));
  %end
  vo.u_cam{ic}{iitt}(:,ill)        = v1D(:,i1D('u_cam'));
  vo.v_cam{ic}{iitt}(:,ill)        = v1D(:,i1D('v_cam'));
  vo.spdt{ic}{iitt}(:,ill)         = v1D(:,i1D('spdt'));
  vo.spdq{ic}{iitt}(:,ill)         = v1D(:,i1D('spdq'));
  vo.camdt{ic}{iitt}(:,ill)        = v1D(:,i1D('camdt'));
  vo.camdq{ic}{iitt}(:,ill)        = v1D(:,i1D('camdq'));
  vo.spmc{ic}{iitt}(:,ill)         = v1D(:,i1D('spmc'));
  vo.spmcup{ic}{iitt}(:,ill)       = v1D(:,i1D('spmcup'));
  %vo.spmcdn{ic}{iitt}(:,ill)       = v1D(:,i1D('spmcdn'));
  vo.thetae{ic}{iitt}(:,ill)       = v1D(:,i1D('thetae'));
  vo.relhum{ic}{iitt}(:,ill)       = v1D(:,i1D('relhum'));
  vo.div{ic}{iitt}(:,ill)          = v1D(:,i1D('div')); 
%  vo.div{ic}{iitt}(:,ill)      =  mydivmax(vo.rho_crm{ic}{iitt}(:,:,ill).*vo.u_crm{ic}{iitt}(:,:,ill), ...
%                                         vo.rho_crm{ic}{iitt}(:,:,ill).*vo.w_crm{ic}{iitt}(:,:,ill),parm.dx,parm.dzi,parm.nzi);
%  vo.div{ic}{iitt}(:,:,ill)    = divergence([vo.u_crm{ic}{iitt}(:,:,ill);vo.u_crm{ic}{iitt}(1,1:nzi,ill)],...
%                                    [vo.w_crm{ic}{iitt}(:,:,ill);vo.w_crm{ic}{iitt}(1,1:nzi,ill)],parm.dx,parm.dzi); % concat the periodic bdry to the other end
%  vo.div{ic}{iitt}(:,ill)       = squeeze(mean(vo.div_crm{ic}{iitt}(mcsllx4Cl{ill},:,ill),1));
  vo.vort{ic}{iitt}(:,:,ill)   = vorticity(vo.u_crm{ic}{iitt}(:,:,ill),vo.w_crm{ic}{iitt}(:,:,ill),parm.dx,parm.dzi);
  vo.enst{ic}{iitt}(ill)       = sum(sum(vo.vort{ic}{iitt}(:,:,ill).^2));
  vo.fftke{ic}{iitt}(:,:,ill)  = myfftke(vo.u_crm{ic}{iitt}(:,:,ill),vo.w_crm{ic}{iitt}(:,:,ill),nx,parm.nzi);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thetae = gettheta_e(p,T,qv)
  Td = tdew(p,qv);
  tlcl = t_lcl(T,Td);
  thetae = theta_e(p,T,qv,tlcl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RH] = specific2relhum(T,P,q)
   qv = q/(1-q);
   qvs = r_sub_s(P,T);
   RH = qv/qvs*100;

function [rho T_rho RH ] = density_temp(T,P,qv,qc,qr,qt,approx)
   % Purpose: use the ideal gas law to get the total density and the
   % associated density tempature for a moist parcel at (T,P)
   Rd   = 287;
   rRd  = 1/Rd;
   Rv   = 461.5; % [J K-1 kg-1]
   RdRv = Rd/Rv; % 0.622
   zvir = Rv/Rd - 1; %Rd/Rv = 461.5/287= 0.622, Rv/Rd=1.606 => zvir=0.61
   % P = 100000 [Pa] 
   %Rd = 287; % [J K-1 kg-1]
   %1/Rd = 0.0034843205
   if (nargout >= 3) % even if T_rho is ~ not called it will be counted as nargout
      qvs = r_sub_s(P,T); % saturation mixing ratio, used for calculating RH = qv/qvs, r_sub_s is used in LI function
      RH  = qv/qvs*100;
   end
   switch approx
   case 1 % for m2005
     T_rho = T.*(1+zvir*qv - qc - qr); % approximate density temperature
   case 2
     T_rho = T*(1+qv/0.622)/(1+qt); % qt = qv+qc+qi
   case 3
     T_rho = T*(1+zvir*qv);
   end
   rho = P.*rRd./T_rho ; % P = rho Ra T_rho => rho = P * 1/Ra * 1/T_rho (1/Rd=0.003484)


function [div divmax_crm] = mydivmax(u,w,dx,dzi,nzi) 
  [ div ] = divergence([u;u(1,1:nzi)],[w;w(1,1:nzi)],dx,dzi); % concat the periodic bdry to the other end
  for iiz = 1:nzi
    [~, imax] = max(abs(div(:,iiz))); % find the maximum divergence at each level over the x-domain, supposing the max represents the MCS convective core
    divmax_crm(iiz) = div(imax,iiz);
  end

function [div] = divergence(u, w, dx, dzi)
  [dudx, ~] = gradient(u,dx,dzi);
  [~, dwdz] = gradient(w,dx,dzi);
  div  = dudx + dwdz;

function [vort] = vorticity(u, w, dx, dzi)
   [~, dudz] = gradient(u,dx,dzi);
   [dwdx, ~] = gradient(w,dx,dzi);
   vort = dudz - dwdx;

function [FFTke]     = myfftke(u, w, nx, nzi)
   FFTu = fft2(u); % fft: columnwise zonal fft for each pressure level independently
   FFTw = fft2(w); % fft2: row and column-wise fft, i.e., fft(fft(uw))
   FFTmagu = abs([FFTu(1,:); 2*FFTu(2:nx/2,:); FFTu(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagw = abs([FFTw(1,:); 2*FFTw(2:nx/2,:); FFTw(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagu = ([FFTmagu(:,1), 2*FFTmagu(:,2:nzi/2), FFTmagu(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTmagw = ([FFTmagw(:,1), 2*FFTmagw(:,2:nzi/2), FFTmagw(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTke = 0.5*(FFTmagu.^2 + FFTmagw.^2);


function [LI, tp500] = liftedIndex(Tp,pp,qv,p,t,tp500)
   % http://www.weathertap.com/guides/aviation/lifted-index-and-k-index-discussion.html
   % http://www.teachingboxes.org/avc/content/Severe_Weather_Indices.htm
   % https://www.ncl.ucar.edu/Support/talk_archives/2010/att-2526/sstats.f
   %  Tp : temp of the parcel
   %  pp : pres of the parcel [Pa]
   %  p  : vertical pressures on sigma-pressure levels [Pa] 
   %  t  : vertical temp on sigma-pressure levels
   %  qv : water vapor mixing ratio of the parcel
   %  lifting index LI = te500 - tp500
   %    Cp     = 1004; % [J kg^-1 K^-1]
   %    Rd     = 287.05; % [J kg^-1 K^-1]
   Rd   = 287;
   Cp   = 1004; %. or 1005.7 % specific heat dry air [J/kg/K]
   RCP  = Rd/Cp; 
   CPR  = Cp/Rd;
   P0   = 100000; % ref pressure for potential temp
   p500 = 50000;
   k500 = findplev(p,p500); % find the closest p-lev to 50000
   te500  = t(k500) + (t(k500+1)-t(k500))/(p(k500+1)-p(k500)) * (p500-p(k500)); % env temp interpolated to closest k500 level, using dry lapse rate
   th     = pottemp(Tp,pp); % potential temp of the parcel 
   Td     = tdew(pp,qv); % dew point temp of the parcel
   tlcl   = t_lcl(Tp,Td); % temp at lcl 
   plcl   = P0 * (tlcl/th)^CPR; % pressure at lcl
   thelcl = theta_e(plcl,tlcl,qv,tlcl);
   tp500  = compT_fr_The(thelcl,p500); % temp of the parcel at p500
   LI  = te500 - tp500; % env - parcel 

function th = pottemp(T,p)
   P0   = 100000; % ref pressure for potential temp
   Rd   = 287;
   Cp   = 1004; %. or 1005.7 % specific heat dry air [J/kg/K]
   RCP  = Rd/Cp; 
   th     = T.*(P0./p).^RCP; % potential temp of the parcel 


function Td = tdew(p,qv)
   % pressure in [Pa], temp in [K]
   RdRv = 0.622; 
   qv = qv+1e-8; 
   e = p*qv/(qv+RdRv);
   loge = log(e); 
   Td = (35.86*loge-4947.2325)/(loge-23.6837);


function Tlcl = t_lcl(Tp,Td)
   % The following code was based on Bolton (1980) eqn #15
   % and claims to have 0.1 K maximum error within -35 < T < 35 C
   %  Tp  = original parcel Temperature in Kelvin
   %  Td  = Temperature at Lifting Condensation Level (K)
   Tlcl = 1.0/(  1.0/(Td-56.0)  + log(Tp/Td)/800 )  + 56.0;


function th_e = theta_e(p, t, qv, tlcl)
   % from widepedia
   % https://en.wikipedia.org/wiki/Equivalent_potential_temperature
   %  e = p*qv/(qv+0.622);
   %  th_l = t*(100000/(p-e))^0.28*(t/tlcl)^(0.28*qv);
   %  th_e = th_l*exp((3036/tlcl-1.78)*qv*(1+0.448*qv));
   % from sstats.f
   qv    = qv + 1e-8;
   power = 0.2854*(1.0 - 0.28*qv);
   xx    = t * (100000.0/p)^power;
   p1    = 3.376/tlcl - 0.00254;
   p2    = (qv*1000.0) * (1.0 + 0.81*qv);
   th_e  = xx*exp(p1*p2);


function tp_p = compT_fr_The(thelcl,p)
   % compute parcel temp at pressure p using theta_e at lcl
   %  p: [Pa]
   %  thelcl: potential temp at LCL 
   %  tp_p: paracel temp at pressure p
   Tguess = (thelcl - 0.5 * max(thelcl-270, 0)^1.05)*(p/100000)^.2;
   epsilon=0.01;
   for iter=1:100
      w1 = r_sub_s(p,Tguess); % saturation mixing ratio
      w2 = r_sub_s(p,Tguess+1);
      tenu = theta_e(p,Tguess,w1,Tguess);
      tenup = theta_e(p,Tguess+1,w2,Tguess+1.);
      cor = (thelcl - tenu) / (tenup - tenu);
      Tguess = Tguess + cor;
      if ( (cor < epsilon) & (-cor < epsilon) ) 
         tp_p = Tguess; return
      end
   end
   thwlcl = theta_wetb(thelcl);
   tp_p = thwlcl*(p/100000.0)^0.286;

function qvs = r_sub_s(p,t)
   % this calls function e_sub_s which computes saturation
   % vapor pressure (Pa) and converts to sat. mixing ratio (kg/kg)
   %  p - pressure (pa)
   %  t  - temperature (k)
   %  qvs : staturation mixing ratio [kg/kg]
   RdRv = 0.622; 
   es = e_sub_s(t);
   qvs = RdRv*es/(p-es); 


function es = e_sub_s(t)
   % compute saturation vapor pressure (Pa) over liquid with
   % polynomial fit of goff-gratch (1946) formulation. (walko, 1991)
   c = [610.5851,44.40316,1.430341,.2641412e-1,.2995057e-3,.2031998e-5,.6936113e-8,.2564861e-11,-.3704404e-13];
   x = max(-80, t-273.16);
   es = c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+x*(c(6)+x*(c(7)+x*(c(8)+x*c(9))))))));


function th_wetb = theta_wetb(thetae)
   % polynomial fit to data in  Smithsonian Meteorological Tables showing Theta-e and Theta-w
   c = [-1.00922292e-10, -1.47945344e-8, -1.7303757e-6, -0.00012709, 1.15849867e-6, -3.518296861e-9, 3.5741522e-12 ];
   d = [0.00000000,   -3.5223513e-10, -5.7250807e-8, -5.83975422e-6, 4.72445163e-8, -1.13402845e-10, 8.729580402e-14];
   x = min(475.0,thetae);
   if ( x <= 335.5 ) 
      th_wetb = 273.15 + c(1)+x*(c(2)+x*(c(3)+x*(c(4)+x*(c(5)+x*(c(6)+ x*c(7) )))));
   else
      th_wetb = 273.15 + d(1)+x*(d(2)+x*(d(3)+x*(d(4)+x*(d(5)+x*(d(6)+ x*d(7) )))));
   end



