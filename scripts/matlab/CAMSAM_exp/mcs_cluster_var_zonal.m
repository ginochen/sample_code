function mcs_cluster_var_zonal
load ~/comp
season = 'JJA';
dlat = 'lat0055';
xpts = 50;
dpts = [num2str(xpts) 'pts'];
dims = '1d';
if strcmp(comp,'PEG')
  casei = 'F_2000_SPCAM_m2005_3hrly1';
  diri = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/F_2000_SPCAM_m2005_3hrly1/atm/hist/';
  fci = [diri casei '.cam.h1'];
  diro = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/F_2000_SPCAM_m2005_3hrly1/atm/hist/' dlat '/' season '/'];
  fmcs = [diro '/mcs_clusters.mat'];
  load([diro '/mcs_cluster_parm.mat']);
  load(fmcs,'t','t4Cl','nt4Cl','mcsn','mcsnll','nCl4nt','nCl','nt4ClRowInd','mcsillt4Cl','mcsllx4Cl','mcsilltcentroids');
  cd([diro '/mcs_cluster_var/' dpts '/' dims])
else
  error('do this on pegasus')
end
if exist('mcs_cluster_zonal_var.mat','file')
  load mcs_cluster_zonal_var
else
  save('mcs_cluster_zonal_var.mat','parm')
end

[Nt,idxc]=sort(cellfun(@length,mcsillt4Cl),'descend'); % obtain the cluster with longest lifecycle in descending order, ic=368 is over US
vo.non = [];
for ic = 51:sum(Nt>3)%4:numel(Nt)
  for it = 1:Nt(ic) % disregard appending the additional 2 time steps at the ends
    fin = [fci '.' t{mcsillt4Cl{idxc(ic)}{it}(1,3)} '.nc'];
    [vi] = readallvar(fin);
    ilat = mcsilltcentroids{idxc(ic)}(it,2);
%      [ilat N] = mode(mcsillt4Cl{idxc(ic)}{it}(:,2)); % ilat: lat with most longitudes, N: # of lon
    ilons = mcsillt4Cl{idxc(ic)}{it}(:,2)==ilat; % ilons associated to the ilat
    ilons = addlons(mcsillt4Cl{idxc(ic)}{it}(ilons,1),parm,xpts); % ilons: add additional lon indices to the detected lons at most repeated latitude
    [vo] = calallvar(it,idxc(ic),ilons,ilat,mcsillt4Cl{idxc(ic)}{it}(1,3),parm,vi,vo);
  end
  [vo] = calallvar_mean_anomaly(idxc(ic),parm,vo);
  disp(ic)
end
save('mcs_cluster_zonal_var','vo','idxc','Nt','-append');



%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vi] = readallvar(fin);
  vi.z3       = flipdim(ncread(fin,'Z3'),3);
  vi.phis     = ncread(fin,'PHIS');
  vi.omega    = flipdim(ncread(fin,'OMEGA'),3);
  vi.u        = flipdim(ncread(fin,'U'),3);
  vi.v        = flipdim(ncread(fin,'V'),3);
  vi.ps       = squeeze(ncread(fin,'PS'));
  vi.T        = flipdim(ncread(fin,'T'),3);
  vi.q        = flipdim(squeeze(ncread(fin,'Q')),3);
  vi.qt       = flipdim(squeeze(ncread(fin,'QT')),3);
  vi.spdt     = flipdim(squeeze(ncread(fin,'SPDT')),3);

function [vo] = calallvar(iitt,ic,ilons,ilat,it,parm,vi,vo)
  for ill=1:numel(ilons)
    l1 = ilons(ill);
    l10 = l1-1; % used for calculating 1st deriv
    l11 = l1+1;
    if l1 == length(parm.lon) % at right bd
      l11 = 1;
    elseif l1 == 1
      l10 = length(parm.lon); % at left bd
    end
    l2 = ilat; % lat
    vo.ilat(iitt) = l2;
    l20 = l2-1;
    l21 = l2+1;
    vo.mcsillt4Cl_fixlat{ic}{iitt}(ill,:) = [l1,l2,it];
  z0 = vi.phis(l1,l2)/parm.g; % surface geopotential height
  z30 = squeeze(vi.z3(l1,l2,1:parm.nz))-z0; % geopotential height w.r.t local surface height
    [pi, pm, dpi, dpm]  = hybrid2p(parm.p0, vi.ps(l1,l2), parm.hyai, parm.hybi, parm.hyam, parm.hybm, parm.nz);
    for iz=1:parm.nz
      rho(iz)  = density_temp(squeeze(vi.T(l1,l2,iz)),pm(iz),vi.q(l1,l2,iz)./(1-vi.q(l1,l2,iz)),[],[],squeeze(vi.qt(l1,l2,iz)), 2); % midpoint density, disregard the rain effect here
    end
  [z, ~] = p2z(rho,dpi,dpm,parm.g,z0,parm.nz);
  disp(z'-z30); 
%  disp(z30); 
pause
    vo.spdt{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)),squeeze(vi.spdt(l1,l2,1:parm.nz))', parm.zint','linear','extrap');
    vo.rho{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)), rho', parm.zint','linear','extrap'); % forgot to use this for pv
    vo.w{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)), -squeeze(vi.omega(l1,l2,1:parm.nz))./(rho'*parm.g), parm.zint','linear','extrap'); % omega to w
    vo.theta{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)),[squeeze(pottemp(squeeze(vi.T(l1,l2,1:parm.nz))',pm))], parm.zint','linear','extrap');
    vo.thetaz{ic}{iitt}(:,ill) = (vo.theta{ic}{iitt}(2:end,ill)-vo.theta{ic}{iitt}(1:end-1,ill))/parm.dzi;
    theta10 = interp1(squeeze(vi.z3(l10,l2,1:parm.nz)),[squeeze(pottemp(squeeze(vi.T(l10,l2,1:parm.nz))',pm))], parm.zint','linear','extrap');
    theta11 = interp1(squeeze(vi.z3(l11,l2,1:parm.nz)),[squeeze(pottemp(squeeze(vi.T(l11,l2,1:parm.nz))',pm))], parm.zint','linear','extrap');
    vo.thetax{ic}{iitt}(:,ill) = (theta11 - theta10)/parm.dlon(l2);
    theta20 = interp1(squeeze(vi.z3(l1,l20,1:parm.nz)),[squeeze(pottemp(squeeze(vi.T(l1,l20,1:parm.nz))',pm))], parm.zint','linear','extrap');
    theta21 = interp1(squeeze(vi.z3(l1,l21,1:parm.nz)),[squeeze(pottemp(squeeze(vi.T(l1,l21,1:parm.nz))',pm))], parm.zint','linear','extrap');
    vo.thetay{ic}{iitt}(:,ill) = (theta21 - theta20)/parm.dlat;
    vo.v_cam{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)),[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear','extrap');
    vo.u_cam{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)),[squeeze(vi.u(l1,l2,1:parm.nz))], parm.zint','linear','extrap');
    vo.uz{ic}{iitt}(:,ill) = (vo.u_cam{ic}{iitt}(2:end,ill)-vo.u_cam{ic}{iitt}(1:end-1,ill))/parm.dzi;
    u20 = interp1(squeeze(vi.z3(l1,l20,1:parm.nz)),[squeeze(vi.u(l1,l20,1:parm.nz))], parm.zint','linear','extrap');
    u21 = interp1(squeeze(vi.z3(l1,l21,1:parm.nz)),[squeeze(vi.u(l1,l21,1:parm.nz))], parm.zint','linear','extrap');
    vo.uy{ic}{iitt}(:,ill) = (u21 - u20)/parm.dlat;
    vo.v_cam{ic}{iitt}(:,ill) = interp1(squeeze(vi.z3(l1,l2,1:parm.nz)),[squeeze(vi.v(l1,l2,1:parm.nz))], parm.zint','linear','extrap');
    vo.vz{ic}{iitt}(:,ill) = (vo.v_cam{ic}{iitt}(2:end,ill)-vo.v_cam{ic}{iitt}(1:end-1,ill))/parm.dzi;
    v10 = interp1(squeeze(vi.z3(l10,l2,1:parm.nz)),[squeeze(vi.v(l10,l2,1:parm.nz))], parm.zint','linear','extrap');
    v11 = interp1(squeeze(vi.z3(l11,l2,1:parm.nz)),[squeeze(vi.v(l11,l2,1:parm.nz))], parm.zint','linear','extrap');
    vo.vx{ic}{iitt}(:,ill) = (v11 - v10)/parm.dlon(l2);
    vo.pv{ic}{iitt}(:,ill) = (-vo.vz{ic}{iitt}(:,ill).*vo.thetax{ic}{iitt}(1:end-1,ill)+vo.uz{ic}{iitt}(:,ill).*vo.thetay{ic}{iitt}(1:end-1,ill)+...
            (parm.co(l2) + vo.vx{ic}{iitt}(1:end-1,ill)-vo.uy{ic}{iitt}(1:end-1,ill)).*vo.thetaz{ic}{iitt}(:,ill))./vo.rho{ic}{iitt}(1:end-1,ill);
    %vo.pv     {ic}{iitt}(:,ill) = vo.pv{ic}{iitt}(:,ill)./vo.rho{ic}{iitt}(1:end-1,ill); % total pv
    %vo.pv0    {ic}{iitt}(:,ill) = parm.co(l2)*vo.thetaz0{ic}./vo.rho0{ic}(1:end-1); % pv MCS mean
  end
  vo.u_storm{ic}{iitt} = mean(mean(vo.u_cam{ic}{iitt}(1:parm.i6km,:),1),2);
  vo.v_storm{ic}{iitt} = mean(mean(vo.v_cam{ic}{iitt}(1:parm.i6km,:),1),2);
  for ill=1:numel(ilons)
    vo.u_stormrel{ic}{iitt}(:,ill) = vo.u_cam{ic}{iitt}(:,ill) - vo.u_storm{ic}{iitt};
    vo.v_stormrel{ic}{iitt}(:,ill) = vo.v_cam{ic}{iitt}(:,ill) - vo.v_storm{ic}{iitt};
  end

function [vo] = calallvar_mean_anomaly(ic,parm,vo)
  vo.rho0{ic}    = mean([vo.rho{ic}{:}],2);
  vo.theta0{ic}  = mean([vo.theta{ic}{:}],2);
  vo.thetaz0{ic} = mean([vo.thetaz{ic}{:}],2);
  for iitt = 1:numel(vo.pv{ic})
    vo.pv0{ic}{iitt}  = parm.co(vo.ilat(iitt))*vo.thetaz0{ic}./vo.rho0{ic}(1:end-1);
    for ill=1:size(vo.pv{ic}{iitt},2)
      vo.pva    {ic}{iitt}(:,ill) = vo.pv{ic}{iitt}(:,ill) - vo.pv0{ic}{iitt}; % pv anomoly from MCS mean
      vo.thetaa{ic}{iitt}(:,ill) = vo.theta{ic}{iitt}(:,ill)-vo.theta0{ic};
    end
  end

function newilons = addlons(ilons,parm,N)
  % extract the lon points around the centroid of MCS
  % N: total # of longitudes to extract (only use even number)
  nlon=numel(parm.lon);
  Nh=ceil(N/2); % half of N
  newilons = zeros(1,N);
  N1=ceil(numel(ilons)/2); % # of 1st half of ilons
  N1c=Nh-N1; % remaining # 1st half of newilons
  N2=numel(ilons)-N1; % # of 2nd half of ilons
  N2c=Nh-N2; % remaining # 2st half of newilons
  newilons(Nh-N1+1:Nh+N2)=ilons; % add ilons to center of the N-element newilons vector
  for i=1:N1c % add to 1st half of newilons
    if ilons(1)-i>0
      newilons(N1c-i+1) = ilons(1)-i;
    else
      newilons(N1c-i+1) = nlon-i+1;
    end
  end
  for i=1:N2c % add to 2nd half of newilons
    if ilons(end)+i<=nlon
      newilons(Nh+N2+i) = ilons(end)+i;
    else
      newilons(Nh+N2+i) = ilons(end)+i-nlon;
    end
  end
