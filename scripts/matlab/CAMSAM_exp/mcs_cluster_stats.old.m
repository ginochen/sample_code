function [ivv] = mcs_cluster_stats(ivv,ilndocn)
%function mcs_cluster_stats(vname,vdim,camcrm{icc})
  % Purpose: Calc mean, median or eof for cluster samples in 1Dspace-time (heating) or 0Dspace-time (precip). 
  %   The input data random variable dimension is the space and time index
  %   (e.g., z-t), and the sample dimension is the clusters. If the space index
  %   is just a point (e.g., surface rain rate), then the r.v. dimension is
  %   just time.
  % Input Arguments:
  %   vname: 'spdt', 'dtcond', 'spdq', 'prect',  etc 
  %   vdim: 0, 1
  %   camcrm: 'CAM', 'CRM'
  %    
% mcs_cluster_stats(1,0)
%ivv = 2 % 2 spdt 3 dtcond 4 spdq 5 prect
%ivv = 1 % precc precl divmax_crm LI_crm FFTke_crm
%ivv = 2 % freq_precc freq_precl du03km dv03km 
%ilndocn = 2 % 0 ('o'cn), 1 ('l'nd), 2 (ocn + lnd), 3 (no condition)
minnt = 2
clustthr = 2;%15; %numel(vtmp2)/nlndocnCl(nt)
%maxnt = 12;
vnames = {'u_crm','w_crm','T_crm','qv_crm','qc_crm','qi_crm','qr_crm','rh_crm','rho_crm','vort','fftke','enst','u_cam','v_cam','spdt','spdq','spmc','spmcup','spmcdn','thetae_crm','relhum_crm',         'mcsillt4Cl','div','LI','precc','precl','freq_precc','freq_precl','du03km','dv03km','ctype','prate'};
  vname = {'crm_var','fprecc_fprecl_du03km_dv03km','spdt','dtcond','spdq','prect'};
display(vname{ivv})
  vdim = [5 6 1 1 1 0]; % 0 (0D) 1 (1D) 5 (crm_var)
  cc   = [2 2 1 1 1 1]; % 1 CAM 2 CRM 
  camcrm = {'CAM','CRM'};
  lndocnstr = {'_ocn','_lnd','_lndocn','nocond'};
  ilndocn = ilndocn + 1; % +1 for actual index
display(lndocnstr{ilndocn})
  spcase      = 'F_2000_SPCAM_m2005_3hrly1';
  spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
  spHist.dir  = [spArchive spcase '/atm/hist/'];
  fdiro = spHist.dir;
  load([spHist.dir 'mcs_clusters'],'t','t4Cl','nt4Cl','mcsnll','nCl4nt','nCl','nt4ClRowInd','mcsillt4Cl');
  if strcmp(lndocnstr{ilndocn},'nocond')
    island = zeros(nCl,1);
    disp('setting island to a zero vector, and make sure to set lndocn = 0 to do all clusters without lnd ocn condition')
  else
    load([spHist.dir 'mcs_clusters'],'island');
  end
  dimstr = {'', ':,', ':,:,'}; % for sprintf to specify 0D 1D and 2D
%  run ~/scripts/matlab/startup.m
  spfile = [spHist.dir spcase '.cam.h1.' t{1} '.nc']; 
  p0   = ncread(spfile,'P0'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  phis = ncread(spfile,'PHIS'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  hyam = flipud(ncread(spfile,'hyam'));
  hybm = flipud(ncread(spfile,'hybm')); 
  hyai = flipud(ncread(spfile,'hyai'));
  hybi = flipud(ncread(spfile,'hybi'));
  nz = 28;
  nx = 32;
  g  = 9.8; 
  fname = 'F_2000_SPCAM_m2005_3hrly1_MCSVarMAM';
  load([spHist.dir fname],'parm')
  if vdim(ivv) == 5
    f1 = sprintf('%s/mcs_stats_CRM_divmax%s.mat',fdiro,lndocnstr{ilndocn});
    f2 = sprintf('%s/mcs_stats_CRM_enst%s.mat',  fdiro,lndocnstr{ilndocn});
    f3 = sprintf('%s/mcs_stats_CRM_FFTke%s.mat', fdiro,lndocnstr{ilndocn});
    f4 = sprintf('%s/mcs_stats_CRM_LI%s.mat',    fdiro,lndocnstr{ilndocn});
    f5 = sprintf('%s/mcs_stats_CRM_precc%s.mat', fdiro,lndocnstr{ilndocn});
    f6 = sprintf('%s/mcs_stats_CRM_precl%s.mat', fdiro,lndocnstr{ilndocn});
    f7 = sprintf('%s/mcs_cluster_extMAM%s.mat',  fdiro,lndocnstr{ilndocn});
    if exist(f1,'file') || exist(f2,'file') || exist(f3,'file') || exist(f4,'file') || exist(f5,'file') || exist(f6,'file') || exist(f7,'file') 
      error('some file exists, exiting')
    else
      disp('saving initial files with only parm attached')
    end
    save(f1,'parm'); save(f2,'parm'); save(f3,'parm'); save(f4,'parm'); save(f5,'parm'); save(f6,'parm'); save(f7,'parm')
  elseif vdim(ivv) == 6
    f8 = sprintf('%s/mcs_stats_CRM_freq_precc%s.mat', fdiro,lndocnstr{ilndocn});
    f9 = sprintf('%s/mcs_stats_CRM_freq_precl%s.mat', fdiro,lndocnstr{ilndocn});
    f10 = sprintf('%s/mcs_stats_CRM_du03km%s.mat',    fdiro,lndocnstr{ilndocn});
    f11 = sprintf('%s/mcs_stats_CRM_dv03km%s.mat',    fdiro,lndocnstr{ilndocn});
    f12 = sprintf('%s/mcs_cluster_extMAM2%s.mat',    fdiro,lndocnstr{ilndocn});
    if exist(f8,'file') || exist(f9,'file') || exist(f10,'file') || exist(f11,'file')  || exist(f12,'file')
      error('some file exists, exiting')
    else
      disp('saving initial files with only parm attached')
    end
    save(f8,'parm'); save(f9,'parm'); save(f10,'parm'); save(f11,'parm'); save(f12,'parm')
  else
    eval(sprintf('save(''%s/mcs_stats_%s_%s%s'',''parm'')',fdiro,camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn},vname{ivv}))
  end
  if ~exist('maxnt','var')
    maxnt = max(nt4Cl);
  end
%  load('F_2000_SPCAM_m2005_3hrly1_MCSVarMAM','varCAM0D_all','varCRM0D_all','varCAM1D_all','varCRM1D_all','varCRM2D_all','stabilityIndex_all','divmax_crm_all','enstrophy_crm_all','FFTke_crm_all','mcs_all')
  % load the interpolated, calculated variables
%  if vdim(ivv)==5 
%    eval(sprintf('load(''%s/%s'',''divmax_crm_all'',''enstrophy_crm_all'',''FFTke_crm_all'',''stabilityIndex_all'',''varCRM0D_all'',''mcs_all'',''parm'',''idv'');',spHist.dir,fname));
%    vtmpdiv = divmax_crm_all; %1D
%    vtmpens = cellfun(@transpose,enstrophy_crm_all,'uniformoutput',false); % transpose all matrix in each cell 
%    vtmpFFTke = FFTke_crm_all;
%    for it = 1:numel(t)
%      vtmpLI{it} = stabilityIndex_all{it}(:,idv.stabilityIndex.LI_crm); % 0D
%      vtmpprecc{it} = varCRM0D_all{it}(:,idv.varCRM0D.precc);
%      vtmpprecl{it} = varCRM0D_all{it}(:,idv.varCRM0D.precl);
%    end
%  elseif vdim(ivv)==6 
%    eval(sprintf('load(''%s/%s'',''varCRM0D_all'',''varCAM0D_all'',''mcs_all'',''parm'',''idv'');',spHist.dir,fname));
%    for it = 1:numel(t)
%      vtmpfprecc{it} = varCRM0D_all{it}(:,idv.varCRM0D.freq_precc);
%      vtmpfprecl{it} = varCRM0D_all{it}(:,idv.varCRM0D.freq_precl);
%      vtmpdu03km{it} = varCRM0D_all{it}(:,idv.varCRM0D.du03km);
%      vtmpdv03km{it} = varCAM0D_all{it}(:,idv.varCAM0D.dv03km);
%    end
%  else
%    eval(sprintf('load(''%s/%s'',''var%s%dD_all'',''mcs_all'',''parm'',''idv'');',spHist.dir,fname,camcrm{cc(ivv)},vdim(ivv)));
%    eval(sprintf('vtmp = var%s%dD_all;',camcrm{cc(ivv)},vdim(ivv)));
%    eval(sprintf('idxv = idv.var%s%dD.%s',camcrm{cc(ivv)},vdim(ivv),vname{ivv}));
%  end
  load([spHist.dir season{iseason} '/mcs_cluster_var/mcs_cluster'],vname{ivv})
  if exist(sprintf('%s/mcs_stats_%s_%s%s.mat',spHist.dir,camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn}),'file'); % check if stats for some nt's are already calc'd, skip them...
    eval(sprintf('load(''%s/%s/mcs_stats_%s_%s%s'',''%s'',''ntvec'',''nlndocnCl'');',spHist.dir,season{iseason},camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn},vname{ivv}));
    ntvec_old = ntvec;
  end
  ntvec = [minnt:maxnt];
  if exist('ntvec_old','var')
    ntvec = setxor(ntvec_old,ntvec); % find the newly added ntvec
    if isempty(ntvec); error('no new nts to do...'); end
    disp(['do these nt''s = ' num2str(ntvec)]);
  end
  divmax=[];enst=[];FFTke=[];LI=[];precc=[];precl=[];freq_precc=[];freq_precl=[]; du03km=[]; dv03km=[];
  nlndocnCl = zeros(numel(ntvec),1);
  for nt = ntvec 
tic
    vtmp2 = []; vtmp3 =[]; vtmp4 = []; vtmp5 = []; vtmp6 = []; vtmp7 = []; vtmp8=[]; vtmp9=[]; vtmp10=[]; vtmp11=[];
    if nCl4nt(nt) % at least one cluster sample exist 
      iiic = 1;
      for iic = 1:nCl4nt(nt)
        ic = nt4ClRowInd{nt}(iic); % ic-cluster with nt = i
        if any(ismember(t4Cl(ic,:),[1 numel(t)])); break; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
        if island(ic)==ilndocn-1  % qualify for ocean index or default no land ocean condition
          for iit = 1:nt
            it = t4Cl(ic,iit); % find  t4Cl(:,1)==it 
            if iit==1  % at the 1st time index save extra -1 index
              spfile = [spHist.dir spcase '.cam.h1.' t{it-1} '.nc']; 
            elseif iit==nt % at the last time index save extra +1 index
              spfile = [spHist.dir spcase '.cam.h1.' t{it+1} '.nc'];  
            end 
            vtmp1=0; irows = 0; precctmp = 0; precltmp = 0; ensttmp = 0; divtmp = 0; FFTketmp = 0; LItmp=0; fprecctmp=0; fprecltmp=0; du03kmtmp=0; dv03kmtmp=0;
            for ill=1:mcsnll(ic,iit) % go through the ilon-ilat at iit time for cluster ic
              irows(ill) = find(ismember(mcs_all{it}.ilonlat, mcsillt4Cl{ic}{iit}(ill,1:2), 'rows')); % find the row of mcs.ilonlat assoc to ill-point at iit-time for ic-cluster
              if ismember(iit,[1,nt])
                ill1 = mcsillt4Cl{ic}{iit}(ill,1:2);
                switch vdim(ivv)
                case 0 
                  vtmp1 = vtmp1 + squeeze(ncread(spfile,upper(vname{ivv}),[ill1(1),ill1(2),1],[1 1 1]));
                case 1
                  z3_cam = squeeze(ncread(spfile,'Z3',[ill1(1),ill1(2),1,1],[1 1 inf 1]));
                  vtmp1  = vtmp1 + interp1(z3_cam, squeeze(ncread(spfile,upper(vname{ivv}),[ill1(1),ill1(2),1,1],[1 1 inf 1])),parm.zint,'linear','extrap');
                case 2
                  vtmp1 = vtmp1 + ncread(spfile,upper(vname{ivv}),[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1]); % (lon, lat, nx, 1, nz)
                case 5 % requires calc variables
tic
                  qc = squeeze(ncread(spfile,'CRM_QC',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
                  prec_crm = squeeze(ncread(spfile,'CRM_PREC',[ill1(1),ill1(2),1,1,1],[1,1,inf,1,1]))*3.6e6;
                  T  = squeeze(ncread(spfile,'T',[ill1(1),ill1(2),1,1],[1 1 inf 1]));
                  Ts = squeeze(ncread(spfile,'TS',[ill1(1),ill1(2),1],[1,1,1]));
                  u  = squeeze(ncread(spfile,'CRM_U',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
                  w  = squeeze(ncread(spfile,'CRM_W',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1]));
                  T_crm  = squeeze(ncread(spfile,'CRM_T',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
                  ps = squeeze(ncread(spfile,'PS',[ill1(1),ill1(2) 1],[1 1 1])); 
                  phis = squeeze(ncread(spfile,'PHIS',[ill1(1),ill1(2) 1],[1 1 1])); 
                  qv = squeeze(ncread(spfile,'CRM_QV',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
toc
                  [pi, pm, dpi, dpm]  = hybrid2p(p0, ps, hyai, hybi, hyam, hybm, nz);
                  for ix = 1:nx
                    for iz = 1:nz
                      [rho(ix,iz), ~, rh_crm(ix,iz)]  = density_temp(T_crm(ix,iz), pm(iz),...
                          qv(ix,iz), 'm2005'); % midpoint density
                    end
                    z0 = phis/g; % surface geopotential height
                    [z, ~] = p2z(rho(ix,:),dpi,dpm,g,z0,nz);
                    uwr(ix,:,:) = interp1(z, [u(ix,:)', w(ix,:)', rho(ix,:)'], parm.zint','linear','extrap');
                  end
                  div0 = mydivmax(uwr(:,:,3).*uwr(:,:,1),uwr(:,:,3).*uwr(:,:,2),parm.dx,parm.dzi,parm.nzi);
                  divtmp = divtmp + div0;
                  vort  = vorticity(uwr(:,:,1),uwr(:,:,2),parm.dx,parm.dzi);
                  enst0 = sum(vort(:).^2);
                  ensttmp = ensttmp + enst0;
                  fftke0 = myfftke(uwr(:,:,1),uwr(:,:,2),nx,parm.nzi);
                  FFTketmp = FFTketmp + fftke0; % use periodogram(ke)
                  LI0 =  liftedIndex(Ts,ps,max(qv(:,1)),pm(2:end),mean(T_crm,1));
                  LItmp = LItmp + LI0;
                  for ix = 1:nx
                    [ctype(ix), prate(ix)] = conv_strat_separation(prec_crm(ix), w(ix,:), qc(ix,:)); 
                  end
                  precc0 = nansum(prate(ctype==1))/nx;
                  precctmp = precctmp + precc0; 
                  precl0 = nansum(prate(ctype==0))/nx;
                  precltmp = precltmp + precl0; 
                  if iit ==1; iitt = 1; elseif iit==nt; iitt = 2;
                  end
                  mcsillt4Cl_ext{ic}{iitt}(ill,1:2) = mcsillt4Cl{ic}{iit}(ill,1:2);
                  div_ext{ic}{iitt}(:,ill) = div0;
                  enst_ext{ic}{iitt}(ill) = enst0;
                  fftke_ext{ic}{iitt}(:,:,ill) = fftke0;
                  LI_ext{ic}{iitt}(ill) = LI0;
                  precc_ext{ic}{iitt}(ill)= precc0;
                  precl_ext{ic}{iitt}(ill)= precl0;
                case 6
%                  z3 = zeros(30,1);
%                  phis = zeros(1,1);
%                  ps = zeros(1,1);
%                  qc = zeros(nx,nz);
%                  qv = zeros(nx,nz);
%                  u = zeros(nx,nz);
%                  w = zeros(nx,nz);
%                  T_crm = zeros(nx,nz);
%                  v = zeros(30,1);
%                  prec_crm = zeros(1,1);
tic
                  z3 = ncread(spfile,'Z3');
                  z3 = squeeze(z3(ill1(1),ill1(2),:));
                  v = ncread(spfile,'V');
                  v = squeeze(v(ill1(1),ill1(2),:));
                  phis = ncread(spfile,'PHIS');
                  phis = squeeze(phis(ill1(1),ill1(2)));
                  qc = ncread(spfile,'CRM_QC');
                  qc = squeeze(qc(ill1(1),ill1(2),:,1,:));
                  qv = ncread(spfile,'CRM_QV');
                  qv = squeeze(qv(ill1(1),ill1(2),:,1,:));
                  prec_crm = ncread(spfile,'CRM_PREC')*3.6e6;
                  prec_crm = squeeze(prec_crm(ill1(1),ill1(2),:,1));
                  u = ncread(spfile,'CRM_U');
                  u = squeeze(u(ill1(1),ill1(2),:,1,:));
                  w = ncread(spfile,'CRM_W');
                  w = squeeze(w(ill1(1),ill1(2),:,1,:));
                  T_crm = ncread(spfile,'CRM_T');
                  T_crm = squeeze(T_crm(ill1(1),ill1(2),:,1,:));
                  ps = ncread(spfile,'PS');
                  ps = squeeze(ps(ill1(1),ill1(2)));
%                  z3 = squeeze(ncread(spfile,'Z3',[ill1(1),ill1(2),1,1],[1 1 inf 1]));
%                  phis = squeeze(ncread(spfile,'PHIS',[ill1(1),ill1(2) 1],[1 1 1])); 
%                  qc = squeeze(ncread(spfile,'CRM_QC',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
%                  prec_crm = squeeze(ncread(spfile,'CRM_PREC',[ill1(1),ill1(2),1,1,1],[1,1,inf,1,1]))*3.6e6;
%                  u  = squeeze(ncread(spfile,'CRM_U',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
%                  w  = squeeze(ncread(spfile,'CRM_W',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1]));
%                  v  = squeeze(ncread(spfile,'V',[ill1(1),ill1(2),1,1],[1 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
%                  T_crm  = squeeze(ncread(spfile,'CRM_T',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
%                  ps = squeeze(ncread(spfile,'PS',[ill1(1),ill1(2) 1],[1 1 1])); 
%                  qv = squeeze(ncread(spfile,'CRM_QV',[ill1(1),ill1(2),1,1,1,1],[1 1 inf 1 inf 1])); % (lon, lat, nx, ny, nz, nt)
toc
                  %close all force
                  v = interp1(z3',[v'], parm.zint','linear','extrap');
                  [pi, pm, dpi, dpm]  = hybrid2p(p0, ps, hyai, hybi, hyam, hybm, nz);
                  for ix = 1:nx
                    for iz = 1:nz
                      [rho(ix,iz), ~, rh_crm(ix,iz)]  = density_temp(T_crm(ix,iz), pm(iz),...
                          qv(ix,iz), 'm2005'); % midpoint density
                    end
                    z0 = phis/g; % surface geopotential height
                    [z, ~] = p2z(rho(ix,:),dpi,dpm,g,z0,nz);
                    uwr(ix,:,:) = interp1(z, [u(ix,:)', w(ix,:)', rho(ix,:)'], parm.zint','linear','extrap');
                  end
                  for ix = 1:nx
                    [ctype(ix), prate(ix)] = conv_strat_separation(prec_crm(ix), w(ix,:), qc(ix,:)); 
                  end
                  fprecc0 = nansum(ctype==1)/nx;
                  fprecctmp = fprecctmp + fprecc0; 
                  fprecl0 = nansum(ctype==0)/nx;
                  fprecltmp = fprecltmp + fprecl0; 
                  du03km0 = diff(mean(uwr(:,[1,parm.i3km],1),1));
                  du03kmtmp = du03kmtmp + du03km0; 
                  dv03km0 = diff(v([1,parm.i3km]));
                  dv03kmtmp = dv03kmtmp + dv03km0; 
                  if iit ==1; iitt = 1; elseif iit==nt; iitt = 2;
                  end
                  mcsillt4Cl_ext{ic}{iitt}(ill,1:2) = mcsillt4Cl{ic}{iit}(ill,1:2);
                  freq_precc_ext{ic}{iitt}(ill) = fprecc0;
                  freq_precl_ext{ic}{iitt}(ill) = fprecl0;
                  du03km_ext{ic}{iitt}(ill) = du03km0;
                  dv03km_ext{ic}{iitt}(ill) = dv03km0;
                  u_ext{ic}{iitt}(:,:,ill)   = uwr(:,:,1);
                  w_ext{ic}{iitt}(:,:,ill)   = uwr(:,:,2);
                  rho_ext{ic}{iitt}(:,:,ill) = uwr(:,:,3);
                  ctype_ext{ic}{iitt}(:,ill) = ctype;
                  prate_ext{ic}{iitt}(:,ill) = prate;
                end  % switch case
              end % if iit == 1 or nt
              %clear functions
            end % ill
            if vdim(ivv) == 5 
              %vtmp2 = myave(vtmp,vtmp1,vtmp2,idxv,irows,nt,it,iit,ic,iiic,vdim,mcsnll) % idxv=1 for all standalone variables
              vtmp2 = myave(vtmpdiv,  divtmp,  vtmp2,1,irows,nt,it,iit,ic,iiic,1,mcsnll);
              vtmp3 = myave(vtmpens,  ensttmp, vtmp3,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp4 = myave(vtmpFFTke,FFTketmp,vtmp4,1,irows,nt,it,iit,ic,iiic,2,mcsnll);
              vtmp5 = myave(vtmpLI,   LItmp,   vtmp5,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp6 = myave(vtmpprecc,precctmp,vtmp6,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp7 = myave(vtmpprecl,precltmp,vtmp7,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
            elseif vdim(ivv) == 6
              vtmp8 = myave(vtmpfprecc,fprecctmp,vtmp8,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp9 = myave(vtmpfprecl,fprecltmp,vtmp9,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp10 = myave(vtmpdu03km,du03kmtmp,vtmp10,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
              vtmp11 = myave(vtmpdv03km,dv03kmtmp,vtmp11,1,irows,nt,it,iit,ic,iiic,0,mcsnll);
            else
              vtmp2 = myave(vtmp,vtmp1,vtmp2,idxv,irows,nt,it,iit,ic,iiic,vdim(ivv),mcsnll);
            end
          end
          iiic = iiic+1; % counter for lnd ocn conditioned cluster
        end
      end   
      nlndocnCl(nt) = iiic-1; % total count of lnd ocn conditioned cluster 
      if vdim(ivv) == 5 
        divmax  = mystats(divmax,vtmp2,nt,clustthr,nlndocnCl);
        enst    = mystats(enst,  vtmp3,nt,clustthr,nlndocnCl);
        FFTke   = mystats(FFTke, vtmp4,nt,clustthr,nlndocnCl);
        LI      = mystats(LI,    vtmp5,nt,clustthr,nlndocnCl);
        precc   = mystats(precc, vtmp6,nt,clustthr,nlndocnCl);
        precl   = mystats(precl, vtmp7,nt,clustthr,nlndocnCl);
        save(f1,'divmax', 'ntvec','nlndocnCl','-append')
        save(f2,'FFTke', 'ntvec','nlndocnCl','-append')
        save(f3,'LI', 'ntvec','nlndocnCl','-append')
        save(f4,'enst', 'ntvec','nlndocnCl','-append')
        save(f5,'precc', 'ntvec','nlndocnCl','-append')
        save(f6,'precl', 'ntvec','nlndocnCl','-append')
        save(f7,'mcsillt4Cl_ext','div_ext','enst_ext','fftke_ext','LI_ext','precc_ext','precl_ext','-append')
      elseif vdim(ivv)==6
        freq_precc   = mystats(freq_precc,vtmp8,nt,clustthr,nlndocnCl);
        freq_precl   = mystats(freq_precl,vtmp9,nt,clustthr,nlndocnCl);
        du03km   = mystats(du03km,vtmp10,nt,clustthr,nlndocnCl);
        dv03km   = mystats(dv03km,vtmp11,nt,clustthr,nlndocnCl);
        save(f8,'freq_precc', 'ntvec','nlndocnCl','-append')
        save(f9,'freq_precl', 'ntvec','nlndocnCl','-append')
        save(f10,'du03km', 'ntvec','nlndocnCl','-append')
        save(f11,'dv03km', 'ntvec','nlndocnCl','-append')
        save(f12,'mcsillt4Cl_ext','freq_precc_ext','freq_precl_ext','du03km_ext','dv03km_ext','u_ext','w_ext','rho_ext','ctype_ext','prate_ext','-append')
      else
        eval(sprintf('%s = mystats(%s,vtmp2,nt,clustthr,nlndocnCl);',vname{ivv},vname{ivv})) 
        eval(sprintf('save(''%s/mcs_stats_%s_%s%s'',''%s'',''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn},vname{ivv}))
      end
    end
toc
  end
  if ismember(vdim(ivv),[5,6]) 
    disp(['saved ' vname{ivv} 's to file its associated mcs_stats_CRM_varname_' lndocnstr{ilndocn} '.mat'])
  else
    disp(['saved ' vname{ivv} ' to file mcs_stats_' camcrm{cc(ivv)} '_' vname{ivv} lndocnstr{ilndocn} '.mat'])
  end
%  if exist('ntvec_old','var') % check if stats for some nt's are already calc'd, skip them...
%    ntvec = sort([ntvec,ntvec_old]); % new and old nt's combined
%    if vdim(ivv) == 5
%      eval(sprintf('save(''%s/mcs_stats_CRM_divmax%s'', ''divmax'', ''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_enst%s'', ''enst'', ''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_FFTke%s'',''FFTke'',''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_LI%s'',   ''LI'',   ''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_precc%s'',''precc'',''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_precl%s'',''precl'',''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,lndocnstr{ilndocn}))
%    else
%      eval(sprintf('save(''%s/mcs_stats_%s_%s%s'',''%s'',''parm'',''ntvec'',''nlndocnCl'',''-append'')',fdiro,camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn},vname{ivv}))
%    end
%  else
%    if vdim(ivv) == 5
%      eval(sprintf('save(''%s/mcs_stats_CRM_divmax%s'', ''divmax'', ''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_enst%s'', ''enst'', ''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_FFTke%s'',''FFTke'',''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_LI%s'',   ''LI'',   ''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_precc%s'',''precc'',''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%      eval(sprintf('save(''%s/mcs_stats_CRM_precl%s'',''precl'',''parm'',''ntvec'',''nlndocnCl'')',fdiro,lndocnstr{ilndocn}))
%    else
%      eval(sprintf('save(''%s/mcs_stats_%s_%s%s'',''%s'',''parm'',''ntvec'',''nlndocnCl'')',fdiro,camcrm{cc(ivv)},vname{ivv},lndocnstr{ilndocn},vname{ivv}))
%    end
%  end





%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [myv] = mystats(myv,vtmp2,nt,clustthr,nlndocnCl)
  % myv: my variable 
  disp('begin median calc')
  myv.median{nt} = median(vtmp2,ndims(vtmp2));
  myv.mean{nt} = mean(vtmp2,ndims(vtmp2));
  disp('end median calc')
  if nlndocnCl(nt) >= clustthr
    disp(['begin pca, mean, variance for time length ' num2str(nt+2) ])
    [myv.svec{nt} myv.variance{nt} myv.pc{nt} myv.pc_mean{nt}] = svdVar(reshape(vtmp2,numel(vtmp2)/nlndocnCl(nt),nlndocnCl(nt))); % varin(nvar,nsamp)
    disp(['finish pca for time length ' num2str(nt+2) ])
  else
    disp('Cluster samples are less than space-time indices, skipping pca...')
  end

function vtmp2 = myave(vtmp,vtmp1,vtmp2,idxv,irows,nt,it,iit,ic,iiic,vdim,mcsnll)
  switch vdim
  case 0
    vtmp2(iit+1,iiic) = mean(vtmp{it}(irows,idxv),vdim+1);
    if iit == 1
      vtmp2(iit,iiic) = vtmp1/mcsnll(ic,iit);
    elseif iit == nt
      vtmp2(iit+2,iiic) = vtmp1/mcsnll(ic,iit);
    end
  case 1
    vtmp2(:,iit+1,iiic) = mean(vtmp{it}(:,irows,idxv),vdim+1);
    if iit == 1
      vtmp2(:,iit,iiic) = vtmp1/mcsnll(ic,iit);
    elseif iit == nt
      vtmp2(:,iit+2,iiic) = vtmp1/mcsnll(ic,iit);
    end
  case 2
    vtmp2(:,:,iit+1,iiic) = mean(vtmp{it}(:,:,irows,idxv),vdim+1);
    if iit == 1
      vtmp2(:,:,iit,iiic) = vtmp1/mcsnll(ic,iit);
    elseif iit == nt
      vtmp2(:,:,iit+2,iiic) = vtmp1/mcsnll(ic,iit);
    end
  end

function [pi pm dpi dpm]  = hybrid2p(p0, ps,hyai,hybi,hyam,hybm,nz)
   % hybrid to p coordinate for CESM outputs
   for iz = 1:nz+1
      pi(iz) = hyai(iz)*p0 + hybi(iz)*ps; % interface pressure
   end
   for iz = 1:nz
      dpi(iz) = pi(iz) - pi(iz+1);
   end
   for iz = 1:nz
      pm(iz) = hyam(iz)*p0 + hybm(iz)*ps; % mid-point pressure
      if (iz > 1)
         dpm(iz) = pm(iz-1) - pm(iz);
      else
         dpm(1) = ps - pm(1);
      end
   end


function [rho T_rho RH ] = density_temp(T,P,qv,compset)
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
   T_rho = T*(1+zvir*qv);
   %T_rho = T*(1+qv/0.622)/(1+qv+qc+qi);
   %T_rho1 = T.*(1+0.61*qv - qc - qi); % approximate density temperature
   rho = P.*rRd./T_rho ; % P = rho Ra T_rho => rho = P * 1/Ra * 1/T_rho (1/Rd=0.003484)

function [z dzi] = p2z(rho,dpi,dpm,g,z0,nz);
   % p to z for CESM outputs
   for iz = 1:nz
      if (iz > 1)
         dzm(iz)  = dpm(iz)/(rho(iz)*g);
         dzi(iz) = dpi(iz)/(rho(iz)*g);
      else
         dzm(1) = dpm(1)/(rho(1)*g) + z0; % dpm = rho * g * dzm from surface to mid-point and so on to the next mid-point, 
                                          % add surface geopotential height phis/g, maybe unecessary since near zero over the ocean 
         dzi(1) = dpi(1)/(rho(1)*g) + z0;
      end
      z(iz) = sum(dzm(1:iz)); % height at each midpoint for CRM
   end;

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


function [divmax_crm] = mydivmax(u,w,dx,dzi,nzi) 
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

function [cloudtype rainrate] = conv_strat_separation(prec, w, qc)
   % cloudtype: conv=1 strat=0
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective-stratiform separation criteria by rainrate  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective criteria (tag the lat lon to see if overlapping events
   % occur) from Braun 2010 paper "Simulation and interpretation of the
   % genesis of Tropical Storm Gert (2005) as part of the NASA tropical
   % cloud systems and processes experiment":
   % if local point precip
   % > 20 mm/hr = 20/1000/3600 m/s = 5.6e-6 m/s <---- rain fall velocity
   % > 2*( averaged 24 neighboring precip rate ), then the 
   % point and its surrounding 24 points are convective (if its just 2D
   % then don't use this criteria 
   rainrate=prec;
   if ( prec > 20 ) % [mm/hr]
      cloudtype = 1; % convective
   else
      % if local point precip hasn't reached surface
      % w > 3 m/s or qliq > 0.5 g/kg = 0.5e-3 kg/kg
      if ( max(w)>3 | max(qc)>0.5e-3 )
         cloudtype = 1;
      elseif ( prec > 0.1 ) %[mm/hr] => 0.1mm/hr=2.4mm/day
         cloudtype = 0;
      else
         cloudtype = NaN;
      end
   end
