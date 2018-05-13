function mcs_cluster_var_xy(ic)
% set variable to vo.var = [], if it exists and need substitution 
% set donew=1 in 1) readallvar()  2) callcalvar()
addpath('/nethome/gchen/scripts/matlab/CAMSAM_exp')
addpath('/nethome/gchen/scripts/matlab/')
donew = 0;
doextrat = 0; % do extra time indices (e.g., before MCS detection, look how PV evolves)
load ~/comp % string used for determining server
season = 'JJA';
dlat = 'lat6060';
%dlat = 'lat2525';
xpts = 30; % use odd number to place the cluster centroid in the center of domain
ypts = xpts;
dpts = [num2str(xpts) 'pts'];
dims = '2d'; 
if strcmp(comp,'PEG')
  casei = 'F_2000_SPCAM_m2005_3hrly2';
  diri = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/' casei '/atm/hist/'];
  fci = [diri casei '.cam.h1'];
  diro = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/' casei '/atm/hist/' dlat '/' season '/'];
  if doextrat
    fo = ['mcs_cluster_xy_var_' sprintf('%02d',ic) '_extrat.mat'];
    % load rho0 theta0 thetaz0 from regular time file (i.e., not extra time), 
    % and use it to calc deviations
    load([diro '/mcs_cluster_var/' dpts '/' dims '/mcs_cluster_xy_var_' sprintf('%02d',ic) '.mat'],'vo')
    tmp.rho0 = vo.rho0; 
    tmp.theta0 = vo.theta0;
    tmp.thetaz0 = vo.thetaz0;
    clear vo
    vo = tmp;
  else
    fo = ['mcs_cluster_xy_var_' sprintf('%02d',ic) '.mat'];
  end
  fmcs = [diro '/mcs_clusters_1.mat'];
  if exist([diro '/mcs_cluster_parm.mat'],'file')
    load([diro '/mcs_cluster_parm.mat']);
  else
    parm = mcs_cluster_parm([fci '.0002-06-01-00000.nc'])
  end
%  load([diro '/mcs_cluster_var/mcs_cluster_u_storm.mat']);
%  load([diro '/mcs_cluster_var/mcs_cluster_v_storm.mat']);
  load(fmcs,'t','t4Cl','nt4Cl','mcsn','mcsnll','nCl4nt','nCl','nt4ClRowInd','mcsillt4Cl','mcsllx4Cl','mcsilltcentroids');
  cd([diro '/mcs_cluster_var/' dpts '/' dims])
else
  error('do this on pegasus')
end
if exist(fo,'file') 
  load(fo)
else
  save(fo,'parm','xpts')
end
if ~exist('vo','var')
  vo = [];
end
%vo.w_precc=[];
%vo.w_precl=[];
%vo.theta0=[];
%vo.thetaz0=[];
[Nt,idxc]=sort(cellfun(@length,mcsillt4Cl),'descend'); % obtain the clusters with longest lifecycle in descending order, ic=368 is over US
%for ic = 1:50%:sum(Nt>3)%4:numel(Nt)
  disp(Nt(ic))
  iic = idxc(ic);
  if doextrat % save extra time 
    itextra=[10:-1:1]; % extra time indices
    for it=1:numel(itextra)
      fin = [fci '.' t{mcsillt4Cl{iic}{1}(1,3)-itextra(it)} '.nc'];
      [vi] = readallvar(fin,donew);
      ilons = addpts(mcsilltcentroids{iic}(1,1),parm,xpts);
      ilats = addpts(mcsilltcentroids{iic}(1,2),parm,ypts);
      [vo] = calallvar(it,iic,ilons,ilats,mcsillt4Cl{iic}{1},parm,vi,vo,donew);
    end
  else 
    for it = 1:Nt(ic) % disregard appending the additional 2 time steps 
tic
%    vo.u_storm{iic}{it} = u_storm{iic}{it+1}; %start at the 2nd index since mcsilltcentroids doesn't include 1st and last additional index 
%    vo.v_storm{iic}{it} = v_storm{iic}{it+1};
      fin = [fci '.' t{mcsillt4Cl{iic}{it}(1,3)} '.nc'];
      [vi] = readallvar(fin,donew);
      ilons = addpts(mcsilltcentroids{iic}(it,1),parm,xpts);
      ilats = addpts(mcsilltcentroids{iic}(it,2),parm,ypts);
%    ilon = mcsilltcentroids{iic}(it,1); % doesn't include 1st and last addition index
%    ilat = mcsilltcentroids{iic}(it,2);
%    ilats = mcsillt4Cl{iic}{it}(:,1)==ilon; % ilats associated to the ilon
%    ilats = addpts(mcsillt4Cl{iic}{it}(ilats,2),parm,ypts); % ilons: add additional lon indices to the detected lons at most repeated latitude
%      [ilat N] = mode(mcsillt4Cl{iic}{it}(:,2)); % ilat: lat with most longitudes, N: # of lon
%    ilons = mcsillt4Cl{iic}{it}(:,2)==ilat; % ilons associated to the ilat
%    ilons = addpts(mcsillt4Cl{iic}{it}(ilons,1),parm,xpts); % ilons: add additional lon indices to the detected lons at most repeated latitude
      [vo] = calallvar(it,iic,ilons,ilats,mcsillt4Cl{iic}{it},parm,vi,vo,donew);
toc
    end
  end
  if ~donew
    [vo] = calallvar_mean_anomaly(iic,parm,vo,doextrat);
  end
  disp(iic)
  save(fo,'vo','idxc','Nt','-append','-v6');
  disp('finished saving')
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vi] = readallvar(fin,donew);
%  vi.z3       = flipdim(ncread(fin,'Z3'),3);
%vi.non = [];
  if donew
    vi.non=[];
    %
  else
    vi.ps       = squeeze(ncread(fin,'PS'));
    vi.T        = flipdim(ncread(fin,'T'),3);
    vi.q        = flipdim(squeeze(ncread(fin,'Q')),3);
    vi.qt       = flipdim(squeeze(ncread(fin,'QT')),3);
    vi.u        = flipdim(ncread(fin,'U'),3);
    vi.v        = flipdim(ncread(fin,'V'),3);
    vi.mpdq    = flipdim(ncread(fin,'MPDQ'), 3);
    vi.macpdq    = flipdim(ncread(fin,'MACPDQ'), 3);
    vi.camdq    = flipdim(ncread(fin,'ZMDQ')+ncread(fin,'EVAPQZM')+ncread(fin,'CMFDQ'),3)+vi.mpdq+vi.macpdq;
    vi.u_crm    = squeeze(ncread(fin,'CRM_U')); 
    vi.w_crm    = squeeze(ncread(fin,'CRM_W')); % surface w_crm = 0, all zeros
    vi.zmdt    = flipdim(ncread(fin,'ZMDT'), 3);
    vi.evaptzm    = flipdim(ncread(fin,'EVAPTZM'), 3);
    vi.cmfdt    = flipdim(ncread(fin,'CMFDT'), 3);
    vi.mpdt    = flipdim(ncread(fin,'MPDT')/1004, 3);
    vi.macpdt    = flipdim(ncread(fin,'MACPDT')/1004, 3);
    vi.zmdq    = flipdim(ncread(fin,'ZMDQ'), 3);
    vi.evapqzm    = flipdim(ncread(fin,'EVAPQZM'), 3);
    vi.cmfdq    = flipdim(ncread(fin,'CMFDQ'), 3);
    vi.qi_crm   = squeeze(ncread(fin,'CRM_QI')); % (lon, lat, nx, ny, nz, nt)
    vi.qc_crm   = squeeze(ncread(fin,'CRM_QC')); % (lon, lat, nx, ny, nz, nt)
    vi.qr_crm   = squeeze(ncread(fin,'CRM_QR')); % (lon, lat, nx, ny, nz, nt)
    vi.qv_crm   = squeeze(ncread(fin,'CRM_QV')); % (lon, lat, nx, ny, nz, nt)
    vi.prec_crm = squeeze(ncread(fin,'CRM_PREC'))*3.6e6; %[mm/hr]
    vi.T_crm    = squeeze(ncread(fin,'CRM_T')); % (lon, lat, nx, ny, nz, nt)
    vi.omega    = flipdim(ncread(fin,'OMEGA'),3);
    vi.Ts       = squeeze(ncread(fin,'TS'));
    vi.camdt    = flipdim(ncread(fin,'ZMDT')+ncread(fin,'EVAPTZM')+ncread(fin,'CMFDT')+(ncread(fin,'MPDT')+ncread(fin,'MACPDT'))/1004, 3);
    vi.spdq     = flipdim(ncread(fin,'SPDQ'),3); % sp mass flux [kg/m^2/s] no difference with DCQ
    vi.spdt     = flipdim(squeeze(ncread(fin,'SPDT')),3);
    vi.cloudtop = flipdim(ncread(fin,'CLOUDTOP'),3);
    vi.concld   = flipdim(ncread(fin,'CONCLD'),3);
    vi.prect     = squeeze(ncread(fin,'PRECT'));
  end

function [vo] = calallvar(iitt,ic,ilons,ilats,mcsillt4Cl,parm,vi,vo,donew)
  it = mcsillt4Cl(1,3);
  xpts = numel(ilons);
  ypts = numel(ilats);
  icent = xpts/2;
  ixp  = icent-8:icent+8;
  iyp  = icent-8:icent+8;
  ixc = xpts/2-14:xpts/2+14;
  iyc = ypts/2-14:ypts/2+14;
  ixcs = xpts/2-5:xpts/2+5;
  iycs = ypts/2-5:ypts/2+5;
  if donew % find centroid direction
    % 7:15 is approx 850hPa to 300hPa
   % vo.ucl{ic}{iitt} = mean(reshape(vi.u(mcsillt4Cl(:,1),mcsillt4Cl(:,2),7:15),[],1)); % cloud layer wind
    %vo.vcl{ic}{iitt} = mean(reshape(vi.v(mcsillt4Cl(:,1),mcsillt4Cl(:,2),7:15),[],1));
    %vo.ujet{ic}{iitt} = mean(reshape(max(vi.u(mcsillt4Cl(:,1),mcsillt4Cl(:,2),1:7),[],3),[],1));
    %vo.vjet{ic}{iitt} = mean(reshape(max(vi.v(mcsillt4Cl(:,1),mcsillt4Cl(:,2),1:7),[],3),[],1));
    %vo.uclust{ic}{iitt} = vo.ucl{ic}{iitt}-vo.ujet{ic}{iitt};
    %vo.vclust{ic}{iitt} = vo.vcl{ic}{iitt}-vo.vjet{ic}{iitt};
  end
  for ilo=1:numel(ilons)
    l1 = ilons(ilo);
    l10 = l1-1; % used for calculating 1st deriv
    l11 = l1+1;
    if l1 == length(parm.lon) % at right bd
      l11 = 1;
    elseif l1 == 1
      l10 = length(parm.lon); % at left bd
    end
    for ila=1:numel(ilats)
      l2 = ilats(ila); % lat
      l20 = l2-1;
      l21 = l2+1;
      if donew
        %z1d = vo.z1d{ic}{iitt}(ilo,ila,:);
        vo.ilons{ic}{iitt} = ilons;
        vo.ilat{ic}{iitt} = ilats;
        [pi, pm, dpi, dpm]  = hybrid2p(parm.p0, vi.ps(l1,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        rho1d  = density_temp(squeeze(vi.T(l1,l2,1:parm.nz)),...
                       pm',squeeze(vi.q(l1,l2,1:parm.nz)./(1-vi.q(l1,l2,1:parm.nz))),...
                       [],[],squeeze(vi.qt(l1,l2,1:parm.nz)), 2); % midpoint density, disregard the rain effect here
        [z1d, dzm1d, dzi1d] = p2z(rho1d,dpi,dpm,parm.g,parm.nz);
        vo.u_cam{ic}{iitt}(ilo,ila,:) = interp1(z1d,[squeeze(vi.u(l1,l2,1:parm.nz))], parm.zint','linear');
        vo.v_cam{ic}{iitt}(ilo,ila,:) = interp1(z1d,[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear');
      else % ~donew
        vo.ilons{ic}{iitt} = ilons;
        vo.ilat{ic}{iitt} = ilats;
        [pi, pm, dpi, dpm]  = hybrid2p(parm.p0, vi.ps(l1,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        rho1d  = density_temp(squeeze(vi.T(l1,l2,1:parm.nz)),...
                                pm',squeeze(vi.q(l1,l2,1:parm.nz)./(1-vi.q(l1,l2,1:parm.nz))),...
                                [],[],squeeze(vi.qt(l1,l2,1:parm.nz)), 2); % midpoint density, disregard the rain effect here
        [z1d, dzm1d, dzi1d] = p2z(rho1d,dpi,dpm,parm.g,parm.nz);
        %
        vo.pm{ic}{iitt}(ilo,ila,:)   = pm; 
        vo.rho1d{ic}{iitt}(ilo,ila,:)   = rho1d; 
        vo.z1d{ic}{iitt}(ilo,ila,:)   = z1d; 
%        z1d = vo.z1d{ic}{iitt}(ilo,ila,:);
        for ip=1:parm.nz
          vo.thetae{ic}{iitt}(ilo,ila,ip) = gettheta_e(vo.pm{ic}{iitt}(ilo,ila,ip),vi.T(l1,l2,ip),vi.q(l1,l2,ip));
        end
        vo.thetaei{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(vo.thetae{ic}{iitt}(ilo,ila,:)),parm.zint','linear');
        for ix=1:parm.nx
          vo.w_crmi{ic}{iitt}(ilo,ila,ix,:)  = interp1(z1d,squeeze(vi.w_crm(l1,l2,ix,1:parm.nz)), parm.zint','linear');
          vo.u_crmi{ic}{iitt}(ilo,ila,ix,:)  = interp1(z1d,squeeze(vi.u_crm(l1,l2,ix,1:parm.nz)), parm.zint','linear');
        end
        vo.zmdt{ic}{iitt}(ilo,ila,:)  = interp1(z1d,squeeze(vi.zmdt(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.evaptzm{ic}{iitt}(ilo,ila,:)=interp1(z1d,squeeze(vi.evaptzm(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.cmfdt{ic}{iitt}(ilo,ila,:)  =interp1(z1d,squeeze(vi.cmfdt(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.mpdt{ic}{iitt}(ilo,ila,:)   =interp1(z1d,squeeze(vi.mpdt(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.macpdt{ic}{iitt}(ilo,ila,:) =interp1(z1d,squeeze(vi.macpdt(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.zmdq{ic}{iitt}(ilo,ila,:)  = interp1(z1d,squeeze(vi.zmdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.evapqzm{ic}{iitt}(ilo,ila,:)=interp1(z1d,squeeze(vi.evapqzm(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.cmfdq{ic}{iitt}(ilo,ila,:)  =interp1(z1d,squeeze(vi.cmfdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.mpdq{ic}{iitt}(ilo,ila,:)   =interp1(z1d,squeeze(vi.mpdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.macpdq{ic}{iitt}(ilo,ila,:) =interp1(z1d,squeeze(vi.macpdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.camdqi{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(vi.camdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.camdq{ic}{iitt}(ilo,ila,:) = vi.camdq(l1,l2,1:parm.nz);
%        vo.w_theta{ic}{iitt}(ilo,ila,:) = interp1(squeeze(vo.theta{ic}{iitt}(ilo,ila,:)),squeeze(vo.w{ic}{iitt}(ilo,ila,:)), [298:370],'linear');
        vo.qvi{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(mean(vi.qv_crm(l1,l2,:,1:parm.nz),3)), parm.zint','linear');
        vo.qii{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(mean(vi.qi_crm(l1,l2,:,1:parm.nz),3)), parm.zint','linear');
        vo.qci{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(mean(vi.qc_crm(l1,l2,:,1:parm.nz),3)), parm.zint','linear');
        vo.qri{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(mean(vi.qr_crm(l1,l2,:,1:parm.nz),3)), parm.zint','linear');
%        vo.w_precc{ic}{iitt}(ilo,ila,:)      = nansum(vi.w_crm(l1,l2,vo.ctype{ic}{iitt}(ilo,ila,:)==1,:),3)/sum(vo.ctype{ic}{iitt}(ilo,ila,:)==1);
%        vo.w_precl{ic}{iitt}(ilo,ila,:)      = nansum(vi.w_crm(l1,l2,vo.ctype{ic}{iitt}(ilo,ila,:)==0,:),3)/sum(vo.ctype{ic}{iitt}(ilo,ila,:)==0);
      
        vo.spdti{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(vi.spdt(l1,l2,1:parm.nz)*86400), parm.zint','linear');
        vo.spdqi{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(vi.spdq(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.camdti{ic}{iitt}(ilo,ila,:) = interp1(z1d,squeeze(vi.camdt(l1,l2,1:parm.nz)*86400), parm.zint','linear');
        vo.spdt{ic}{iitt}(ilo,ila,:) = vi.spdt(l1,l2,1:parm.nz)*86400;
        vo.spdq{ic}{iitt}(ilo,ila,:) = vi.spdq(l1,l2,1:parm.nz);
        vo.camdt{ic}{iitt}(ilo,ila,:) = vi.camdt(l1,l2,1:parm.nz)*86400;
        rhos  = density_temp(squeeze(vi.Ts(l1,l2)),...
                                vi.ps(l1,l2),squeeze(vi.q(l1,l2,1)./(1-vi.q(l1,l2,1))),...
                                [],[],squeeze(vi.qt(l1,l2,1)), 2); % midpoint density, disregard the rain effect here
        vo.ws{ic}{iitt}(ilo,ila) = -squeeze(vi.omega(l1,l2,1))./(rhos'*parm.g); % omega to w
        vo.thetas{ic}{iitt}(ilo,ila) = squeeze(pottemp(squeeze(vi.Ts(l1,l2))',vi.ps(l1,l2)));
        vo.Ts{ic}{iitt}(ilo,ila) = vi.Ts(l1,l2);
        [~,idx]=max(vi.cloudtop(l1,l2,1:parm.nz),[],3);
        vo.ctoph{ic}{iitt}(ilo,ila) = z1d(idx); 
        icc = find(squeeze(vi.concld(l1,l2,1:parm.nz))~=0);
        if ~isempty(icc)
          vo.cthick{ic}{iitt}(ilo,ila) = z1d(max(icc))-z1d(min(icc)); % vert thickness of cloud
        else
          vo.cthick{ic}{iitt}(ilo,ila) = 0;
        end 
        [ctype, prate] = conv_strat_separation(squeeze(vi.prec_crm(l1,l2,:)), squeeze(vi.w_crm(l1,l2,:,:)), squeeze(vi.qc_crm(l1,l2,:,:)), squeeze(vi.T_crm(l1,l2,:,:)),pm); 
%        [ctype, prate] = conv_strat_separation(squeeze(vi.prec_crm(l1,l2,:)), squeeze(vi.w_crm(l1,l2,:,:)), squeeze(vi.qc_crm(l1,l2,:,:)), squeeze(vi.T_crm(l1,l2,:,:)),vo.pm{ic}{iitt}(ilo,ila,:));  % use vo.pm instead of pm to speed up
%        vo.ctype{ic}{iitt}(ilo,ila,:)    = ctype;
        vo.frac_precl{ic}{iitt}(ilo,ila) = nansum(prate(ctype==0))/nansum(prate);
        vo.frac_prec{ic}{iitt}(ilo,ila,1)= nansum(prate(prate>=0.1))/nansum(prate);
        vo.frac_prec{ic}{iitt}(ilo,ila,2)= nansum(prate(prate>=3))/nansum(prate);
        vo.frac_prec{ic}{iitt}(ilo,ila,3)= nansum(prate(prate>=5))/nansum(prate);
        vo.frac_precc{ic}{iitt}(ilo,ila) = nansum(prate(ctype==1))/nansum(prate);
        vo.freq_precc{ic}{iitt}(ilo,ila) = nansum(ctype==1)/parm.nx;
        vo.freq_precl{ic}{iitt}(ilo,ila) = nansum(ctype==0)/parm.nx;
        vo.prate{ic}{iitt}(ilo,ila,:)    = prate;
        vo.precc{ic}{iitt}(ilo,ila)      = nansum(prate(ctype==1))/sum(ctype==1);
        vo.precl{ic}{iitt}(ilo,ila)      = nansum(prate(ctype==0))/sum(ctype==0);
        vo.u{ic}{iitt}(ilo,ila,:) = squeeze(vi.u(l1,l2,1:parm.nz));
        vo.v{ic}{iitt}(ilo,ila,:) = squeeze(vi.v(l1,l2,1:parm.nz));
        u_cam = interp1(z1d, squeeze(vi.u(l1,l2,1:parm.nz)), parm.zint','linear');
        v_cam = interp1(z1d, squeeze(vi.v(l1,l2,1:parm.nz)), parm.zint','linear');
        vo.dU03km_cam{ic}{iitt}(ilo,ila) = sqrt(diff(u_cam([1 parm.i3km]))^2+diff(v_cam([1 parm.i3km]))^2);
        vo.du{ic}{iitt}(ilo,ila) = diff(u_cam([1 parm.i3km]));
        [rho2d, T_rho]  = density_temp(squeeze(vi.T_crm(l1,l2,:,1:parm.nz)), ...
                          repmat(pm(1,1:parm.nz),[parm.nx,1]),...
                          squeeze(vi.qv_crm(l1,l2,:,1:parm.nz)),...
                          squeeze(vi.qc_crm(l1,l2,:,1:parm.nz)),...
                          squeeze(vi.qr_crm(l1,l2,:,1:parm.nz)),[], 1); % midpoint density
        th_rho = pottemp(T_rho,repmat(pm(1,1:parm.nz),[parm.nx,1]));
        th_rho_ave=mean(th_rho,1); % LS ave
        for ix = 1:parm.nx
          [z2d(ix,:), dzm2d(ix,:), dzi2d(ix,:)] = p2z(rho2d(ix,:),dpi,dpm,parm.g,parm.nz);
          [B(ix) buoy(ix,:)] = cold_pool_intensity(z2d(ix,:)',th_rho(ix,:),th_rho_ave,parm);
          LI(ix) = liftedIndex(vi.Ts(l1,l2),vi.ps(l1,l2),vi.qv_crm(l1,l2,ix,1),pm,squeeze(vi.T(l1,l2,:)));
        end
%        vo.buoy{ic}{iitt}(:,:,ill)   = buoy; % buoyancy 2D profile
%        vo.B{ic}{iitt}(ilo,ila,:)        = B;
        vo.Bave{ic}{iitt}(ilo,ila)       = mean(B); % buoyancy strength, originally vo.B{ic}{iitt}(:,ill) <- not much use saving this
        vo.Bmax{ic}{iitt}(ilo,ila)       = max(B); % buoyancy strength, originally vo.B{ic}{iitt}(:,ill) <- not much use saving this
        vo.LIave{ic}{iitt}(ilo,ila)      = mean(LI);
        vo.LImax{ic}{iitt}(ilo,ila)      = max(LI);
        vo.LImin{ic}{iitt}(ilo,ila)      = min(LI);
%        vo.LI{ic}{iitt}(ilo,ila,:)       = LI;
        if ilo==1 & ila==1
          vo.u_storm{ic}{iitt} = mean(mean(mean(vi.u(ilons(ixc),ilats(iyc),1:12)))); % 12 is close to the 6km p-level
          vo.v_storm{ic}{iitt} = mean(mean(mean(vi.v(ilons(ixc),ilats(iyc),1:12)))); % 12 is close to the 6km p-level
        end
        vo.mcsillt4Cl_xy{ic}{iitt}(ilo,ila,:) = [l1,l2,it];
        [pi0, pm0, dpi0, dpm0]  = hybrid2p(parm.p0, vi.ps(l1,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        [pi10, pm10, dpi10, dpm10]  = hybrid2p(parm.p0, vi.ps(l10,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        [pi11, pm11, dpi11, dpm11]  = hybrid2p(parm.p0, vi.ps(l11,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        [pi20, pm20, dpi20, dpm20]  = hybrid2p(parm.p0, vi.ps(l1,l20), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        [pi21, pm21, dpi21, dpm21]  = hybrid2p(parm.p0, vi.ps(l1,l21), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
        vo.ps{ic}{iitt}(ilo,ila) = vi.ps(l1,l2);
        vo.prect{ic}{iitt}(ilo,ila) = vi.prect(l1,l2);
        for iz=1:parm.nz
          rho0(iz)  = density_temp(squeeze(vi.T(l1,l2,iz)),pm0(iz),vi.q(l1,l2,iz)./(1-vi.q(l1,l2,iz)),[],[],squeeze(vi.qt(l1,l2,iz)), 2); % midpoint density, disregard the rain effect here
          rho10(iz)  = density_temp(squeeze(vi.T(l10,l2,iz)),pm10(iz),vi.q(l10,l2,iz)./(1-vi.q(l10,l2,iz)),[],[],squeeze(vi.qt(l10,l2,iz)), 2);
          rho11(iz)  = density_temp(squeeze(vi.T(l11,l2,iz)),pm11(iz),vi.q(l11,l2,iz)./(1-vi.q(l11,l2,iz)),[],[],squeeze(vi.qt(l11,l2,iz)), 2);
          rho20(iz)  = density_temp(squeeze(vi.T(l1,l20,iz)),pm20(iz),vi.q(l1,l20,iz)./(1-vi.q(l1,l20,iz)),[],[],squeeze(vi.qt(l1,l20,iz)), 2);
          rho21(iz)  = density_temp(squeeze(vi.T(l1,l21,iz)),pm21(iz),vi.q(l1,l21,iz)./(1-vi.q(l1,l21,iz)),[],[],squeeze(vi.qt(l1,l21,iz)), 2);
        end
        [z0, dzm1d, dzi1d] = p2z(rho0,dpi0,dpm0,parm.g,parm.nz);
        [z10, dzm1d, dzi1d] = p2z(rho10,dpi10,dpm10,parm.g,parm.nz);
        [z11, dzm1d, dzi1d] = p2z(rho11,dpi11,dpm11,parm.g,parm.nz);
        [z20, dzm1d, dzi1d] = p2z(rho20,dpi20,dpm20,parm.g,parm.nz);
        [z21, dzm1d, dzi1d] = p2z(rho21,dpi21,dpm21,parm.g,parm.nz);
%        vo.spdt{ic}{iitt}(ilo,ila,:) = interp1(z0,squeeze(vi.spdt(l1,l2,1:parm.nz)*86400)', parm.zint','linear');
%        plot(parm.zint',squeeze(vo.spdt{ic}{iitt}(ilo,ila,:)),'*'); hold on
%        plot(z0,squeeze(vi.spdt(l1,l2,1:parm.nz))*86400,'r'); hold off
%paauuse
        vo.rho{ic}{iitt}(ilo,ila,:) = interp1(z0, rho0', parm.zint','linear'); % forgot to use this for pv
        vo.w{ic}{iitt}(ilo,ila,:) = interp1(z0, -squeeze(vi.omega(l1,l2,1:parm.nz))./(rho0'*parm.g), parm.zint','linear'); % omega to w
        vo.theta{ic}{iitt}(ilo,ila,:) = interp1(z0,[squeeze(pottemp(squeeze(vi.T(l1,l2,1:parm.nz))',pm0))], parm.zint','linear');
        vo.thetaz{ic}{iitt}(ilo,ila,:) = (vo.theta{ic}{iitt}(ilo,ila,2:end)-vo.theta{ic}{iitt}(ilo,ila,1:end-1))/parm.dzi;
        theta10 = interp1(z10,[squeeze(pottemp(squeeze(vi.T(l10,l2,1:parm.nz))',pm10))], parm.zint','linear');
        theta11 = interp1(z11,[squeeze(pottemp(squeeze(vi.T(l11,l2,1:parm.nz))',pm11))], parm.zint','linear');
        vo.thetax{ic}{iitt}(ilo,ila,:) = (theta11 - theta10)/parm.dlon(l2); % dlon is the lon grid point width in meters at different lats
        theta20 = interp1(z20,[squeeze(pottemp(squeeze(vi.T(l1,l20,1:parm.nz))',pm20))], parm.zint','linear');
        theta21 = interp1(z21,[squeeze(pottemp(squeeze(vi.T(l1,l21,1:parm.nz))',pm21))], parm.zint','linear');
        vo.thetay{ic}{iitt}(ilo,ila,:) = (theta21 - theta20)/parm.dlat;
        vo.v_cam{ic}{iitt}(ilo,ila,:) = interp1(z0,[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear');
        vo.u_cam{ic}{iitt}(ilo,ila,:) = interp1(z0,[squeeze(vi.u(l1,l2,1:parm.nz))], parm.zint','linear');
        vo.uz{ic}{iitt}(ilo,ila,:) = (vo.u_cam{ic}{iitt}(ilo,ila,2:end)-vo.u_cam{ic}{iitt}(ilo,ila,1:end-1))/parm.dzi;
        u20 = interp1(z20,[squeeze(vi.u(l1,l20,1:parm.nz))], parm.zint','linear');
        u21 = interp1(z21,[squeeze(vi.u(l1,l21,1:parm.nz))], parm.zint','linear');
        vo.uy{ic}{iitt}(ilo,ila,:) = (u21 - u20)/parm.dlat;
%        vo.v_cam{ic}{iitt}(ilo,ila,:) = interp1(z0,[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear');
        vo.vz{ic}{iitt}(ilo,ila,:) = (vo.v_cam{ic}{iitt}(ilo,ila,2:end)-vo.v_cam{ic}{iitt}(ilo,ila,1:end-1))/parm.dzi;
        v10 = interp1(z10,[squeeze(vi.v(l10,l2,1:parm.nz))], parm.zint','linear');
        v11 = interp1(z11,[squeeze(vi.v(l11,l2,1:parm.nz))], parm.zint','linear');
        vo.vx{ic}{iitt}(ilo,ila,:) = (v11 - v10)/parm.dlon(l2);
        vo.pv{ic}{iitt}(ilo,ila,:) = (-vo.vz{ic}{iitt}(ilo,ila,:).*vo.thetax{ic}{iitt}(ilo,ila,1:end-1)+vo.uz{ic}{iitt}(ilo,ila,:).*vo.thetay{ic}{iitt}(ilo,ila,1:end-1)+...
                (parm.co(l2) + vo.vx{ic}{iitt}(ilo,ila,1:end-1)-vo.uy{ic}{iitt}(ilo,ila,1:end-1)).*vo.thetaz{ic}{iitt}(ilo,ila,:))./vo.rho{ic}{iitt}(ilo,ila,1:end-1);
        %vo.pv     {ic}{iitt}(ilo,ila,:) = vo.pv{ic}{iitt}(ilo,ila,:)./vo.rho{ic}{iitt}(1:end-1); % total pv
        %vo.pv0    {ic}{iitt}(ilo,ila,:) = parm.co(l2)*vo.thetaz0{ic}./vo.rho0{ic}(1:end-1); % pv MCS mean
        vo.u_stormrel{ic}{iitt}(ilo,ila,:) = vo.u_cam{ic}{iitt}(ilo,ila,:) - vo.u_storm{ic}{iitt};
        vo.v_stormrel{ic}{iitt}(ilo,ila,:) = vo.v_cam{ic}{iitt}(ilo,ila,:) - vo.v_storm{ic}{iitt};
      end
    end
  end
  vo = envwind(vo,ic,iitt,ixp,iyp,parm);


function [vo] = calallvar_mean_anomaly(ic,parm,vo,doextrat)
  if ~doextrat % if doextrat then the rho0, theta0, thetaz0 are all loaded already
    vo.rho0{ic}    = squeeze(mean(mean([vo.rho{ic}{:}],1),2)); % xy env mean over lifetime of MCS
    vo.theta0{ic}  = squeeze(mean(mean([vo.theta{ic}{:}],1),2));
    vo.thetaz0{ic} = squeeze(mean(mean([vo.thetaz{ic}{:}],1),2));
  end
%  vo.rho0{ic}{iitt}    = squeeze(mean(mean([vo.rho{ic}{iitt}],1),2)); % storm xy mean at fixed it
%  vo.theta0{ic}{iitt}  = squeeze(mean(mean([vo.theta{ic}{iitt}],1),2));
%  vo.thetaz0{ic}{iitt} = squeeze(mean(mean([vo.thetaz{ic}{iitt}],1),2));
  for iitt = 1:numel(vo.pv{ic})
    for ila=1:size(vo.pv{ic}{iitt},2)
      vo.pv0{ic}{iitt}(ila,:) = parm.co(vo.ilat{ic}{iitt}(ila))*vo.thetaz0{ic}./vo.rho0{ic}(1:end-1);
      for ilo=1:size(vo.pv{ic}{iitt},1)
        vo.pva{ic}{iitt}(ilo,ila,:) = squeeze(vo.pv{ic}{iitt}(ilo,ila,:)) - vo.pv0{ic}{iitt}(ila,:)'; % pv anomoly from MCS mean
        vo.thetaa{ic}{iitt}(ilo,ila,:) = squeeze(vo.theta{ic}{iitt}(ilo,ila,:)) - vo.theta0{ic};
      end
    end
  end

function newipts = addpts(ilons,parm,N)
  % Add new lon around the existing MCS lon to fill the total N lon 
  % N: total # of longitudes to extract (only use even number)
  nlon=numel(parm.lon);
  Nh=ceil(N/2); % half of N
  newipts = zeros(1,N);
  N1=ceil(numel(ilons)/2); % # of 1st half of ilons
  N1c=Nh-N1; % remaining # 1st half of newilons
  N2=numel(ilons)-N1; % # of 2nd half of ilons
  N2c=Nh-N2; % remaining # 2st half of newilons
  newipts(Nh-N1+1:Nh+N2)=ilons; % add ilons to center of the N-element newilons vector
  for i=1:N1c % add to 1st half of newilons
    if ilons(1)-i>0
      newipts(N1c-i+1) = ilons(1)-i;
    else
      newipts(N1c-i+1) = nlon-i+1;
    end
  end
  for i=1:N2c % add to 2nd half of newilons
    if ilons(end)+i<=nlon
      newipts(Nh+N2+i) = ilons(end)+i;
    else
      newipts(Nh+N2+i) = ilons(end)+i-nlon;
    end
  end


function vo = envwind(vo,ic,iit,ixp,iyp,parm)
 % calc the vert env wind profile to plot on the right y-axis in vertcross figures
 % up: projected u,v vector onto the mcs propag direction vector
 % uenv: the vert env wind profile
  lonc = parm.lon(vo.ilons{ic}{iit});
  latc = parm.lat(vo.ilat{ic}{iit});
  umcs = mean(mean(mean(vo.u_cam{ic}{iit}(ixp,iyp,1:parm.i6km))));% synoptic low propagation speed
  vmcs = mean(mean(mean(vo.v_cam{ic}{iit}(ixp,iyp,1:parm.i6km))));
  vo.umcs{ic}{iit} = umcs; 
  vo.vmcs{ic}{iit} = vmcs; 
  pvecmcs=-[umcs,vmcs]/norm([umcs,vmcs]); % proj storm unit vector, neg for inflow
  [ilonlat,xsquall]=findsquall(lonc,latc,-pvecmcs');
  vo.ilonlat{ic}{iit} = ilonlat;
  vo.xsquall{ic}{iit} = xsquall;
  for ij=1:size(ilonlat,1)
    uzmcs(ij,:)      = vo.u_cam{ic}{iit}(ilonlat(ij,1),ilonlat(ij,2),:)-umcs; % horiz mcs rel flow
    vzmcs(ij,:)      = vo.v_cam{ic}{iit}(ilonlat(ij,1),ilonlat(ij,2),:)-vmcs;  
    up(ij,:)      = [uzmcs(ij,:)',vzmcs(ij,:)']*pvecmcs'; % project the u,v vector onto the unit vector
  end
  vo.uenv{ic}{iit} = mean(up,1);
  vo.up{ic}{iit} = up; % horiz and vert section over the mcs prog dir-plane
  vo.uzmcs{ic}{iit} = uzmcs;
  vo.vzmcs{ic}{iit} = vzmcs;
