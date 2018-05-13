function getBigQPCVar(it,fname)
% function getBigQPCVar(it,fname)
% Find the greatest mismatching heating points for CAM vs SPCAM
%
% - (originally called getBigQPC1Var, ke_spectrum, then divVort, now getBigQPCVar)
% - extract and vertically interpolate (p-to-km) the CRM variables from
%   SPCAM on the large PC1 vertical heating profile space-time indices
%
% save interpolated variables in the order 'u' 'w' 'T' '
% 
% MCS clustering method: 
%    clusterdata(): 
%        heirarchical clustering the convective 
%        gridpoints into MCS clusters
%        uses pdist='euclidean',linkagemethod='single',
%        'criterion'='distance=2.8degrees_lat',cutoff
%    kmedoids(): 
%        search for the center lon-lat for each cluster uses
%        'algorithm'='small', pdist='euclidean'
%
% Save the CAM stratiform precl (large-scale prec rate) on ZM-deep_conv points
run /nethome/gchen/scripts/matlab/startup.m
fname
%%%%%%%%%%%%%%%%%%%%%%%%
%    Job conditions    %
%%%%%%%%%%%%%%%%%%%%%%%%
dofft=0;
dodivvort=0;
docape=0;
%%%%%%%%%%%%%%%%%%%%%
%      Paths        %
%%%%%%%%%%%%%%%%%%%%%
% CAUTION!!! 
% 'it=1' starts at 03600 instead of 01800,
% due to *.rh0.* needed minus 1 from previous timestep.
% spcam_actual_m2005_f09f09_branch.cam.rh0.*.nc and
% F_2000_4SPCAM_m200501.cam.rh0.*.nc both start from 0001-01-14-03600 unlike 
% spcam_actual_m2005_f09f09_branch.cam.r.*.nc starts from 0001-01-14-01800 
compset           = 'm2005';
spArchive         = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
archive           = ['/projects/rsmas/kirtman/gchen/archive/'];
spcamCaser.dir    = [spArchive 'spcam_actual_m2005_f09f09_branch/atm/rest/'];
%Case.dir    = '/projects/rsmas/kirtman/gchen/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
spcamdir_r        = dir([spcamCaser.dir '/spcam_actual_m2005_f09f09_branch.cam.r.*-*-*-*.nc']);
spcamCaser.name   = {spcamdir_r.name};
camCaser.dir      = [archive 'F_2000_4SPCAM_m200501/atm/rest/']; 
camCaserh0.dir    = camCaser.dir;
camCaser.name     = 'F_2000_4SPCAM_m200501.cam.r';
camCaserh0.name   = [camCaser.name 'h0'];
spcamCaserh0.dir  = [archive 'spcam_cam_rh0_m2005/atm/rest/']; % minus previous timestep
spcamCaserh0.name = 'spcam_actual_m2005_f09f09_branch.cam.rh0';
spcamCaseh1.dir  = [spArchive 'spcam_actual_m2005_f09f09_branch_4CRMPREC/atm/hist/']; % minus previous timestep
spcamCaseh1.name = 'spcam_actual_m2005_f09f09_branch_4CRMPREC.cam.h1'; % contains 'CRM_PREC', 'CRM_QRL', 'CRM_QRS'
landfrac          = loadvar('LANDFRAC',camCaserh0,1,0);
noLand = 0; % condition to discard land points, 1 (discard land points), 0 (include land points)
%%%%%%%%%%%%%%%%%%%%%%%%
%    Load dimensions   %
%%%%%%%%%%%%%%%%%%%%%%%%
CaseDim.name   = 'spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600'; % make sure to start at 03600 instead of 01800, and use only this file for the dimensions 
CaseDim.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
p0   = loaddim('P0',CaseDim); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
phis = loaddim('PHIS',CaseDim); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
hyam = flipud(loaddim('hyam',CaseDim));
hybm = flipud(loaddim('hybm',CaseDim)); 
hyai = flipud(loaddim('hyai',CaseDim));
hybi = flipud(loaddim('hybi',CaseDim));
lat  = loaddim('lat',CaseDim);
lon  = loaddim('lon',CaseDim);
lev  = loaddim('lev',CaseDim);  % level centers
ilev = loaddim('ilev',CaseDim); % level interfaces, two interfaces surround one level center
                                % index 1 (TOA pressure) to nlev (surface pressure) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Global Const (global to this function only)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rd   = 287;
rRd  = 1/Rd;
Rv   = 461.5; % [J K-1 kg-1]
RdRv = Rd/Rv; % 0.622
zvir = Rv/Rd - 1; %Rd/Rv = 461.5/287= 0.622, Rv/Rd=1.606 => zvir=0.61
Cp   = 1004; %. or 1005.7 % specific heat dry air [J/kg/K]
RCP  = Rd/Cp; 
CPR  = Cp/Rd;
g    = 9.8; 
P0   = 100000; % ref pressure for potential temp
p500 = 50000; % 50000 Pa
Lv   = 2.5104e6; %; [J/kg]=[m2/s2]  Latent Heat of Vaporization of Water
%%%%%%%%%%%%%%%%%%
%   parameters   %
%%%%%%%%%%%%%%%%%%
nlat = numel(lat); %192
nlon = numel(lon); %288
nlev = numel(lev);
nx    = 32; % evenly spaced
nz    = 28; % levels of CAM from surface to nz levels (according to cam/crm_physics.F90 line 992)
km2m  = 1000; % 1km = 1000m
dx    = 4*km2m; % 4 km grid box
nlons_half = 5; % half the number of indices surrounding the big PC1 longitude, total is nlons_half*2
nlons = nlons_half*2+1;
%nlag  = 7; % extract the nlag timesteps for each grid point that exceeds the PC1 threshold at 0-lag timestep
nlag = 40;
lagt = nlag:-10:-nlag;%nlag:-1:-nlag; % actual lag index 
ilag0  = find(lagt==0); % 0lag index
ilagn0 = find(lagt~=0); % non-0lag index
lagt = [lagt(ilag0) lagt(ilagn0)];% put 0-lag at the first index to evaluate stability and MCS centroid
ilagt = [ilag0,ilagn0]; % 
nlagt = numel(lagt); 
%lagt_extra = nlag+1:-1:-nlag; % extra lag is to get the instantaneous field by stracting previous timestep
% dtlag = 30; % 1800 sec = 30 min for CAM timestep
% timelag = nlag*dtlag:-dtlag:0; % 0 lag is the big PC1 event timestep
nzi  = 10; % # of interpolated z levels, nzi=61(~500m), nzi=180(~150m)
maxz = 18000;
minz = 250; 
zint = linspace(minz,maxz,nzi)'; % interpolated height with ~500m grid-spacing
dzi  = zint(2)-zint(1);
i3km = find(abs(zint-3000) == min(abs(zint-3000)));
dz3km   = zint(i3km) - zint(1); % the distance between first level to 3km 
iLI=[1]; % specify the parcel-levels for lifted index calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Heating profile SVD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMexcUW_UWexcZM_basinwise-svd.mat'
threshold_bigPC = 10; % the 1st EOF mode largest magnitude vertically is 0.5, 
                      % and the average PC has a vertical maximum of 5 K/day heating rate, 
                      % so a threshold of 10 to get 10*0.5=5 is necessary
iEOF = 1; % EOF1 index
%start = [1 1 1]; % starting index for lon(1:288) lat(1:192) nx*nz(1:896) (maybe use lat(76:117) 20 degrees north-south?) 
%count = [inf inf inf]; % number of index to read in var, inf means use all elements in that dimension
bsign = [1 -1 1 1]; % the sign for PC1 in different basin since basin 2 has a different signed distribution
for itt=1:nlagt % lagging big PC1's index by ilag, so ilag=0 is the big PC1 time index, negative values are lead
tic
   ilag=lagt(itt);
   itlag=ilagt(itt);
   %%%%%%%%%%%%%%%%%%%%%%%
   % Load CRM Variables  %
   %%%%%%%%%%%%%%%%%%%%%%%
   ps = ncread([spcamCaser.dir spcamCaser.name{it-ilag}],'PS',[1 1], [inf inf]); %[Pa] in SI units not [hPa] 
   u_crm  = load_reshape(spcamCaser, it-ilag, 'CRM_U',  [1 1 1], [inf inf inf], nlon, nlat, nx, nz); % nz=28
   w_crm  = load_reshape(spcamCaser, it-ilag, 'CRM_W',  [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   T_crm  = load_reshape(spcamCaser, it-ilag, 'CRM_T',  [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qT_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QT', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qc_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QC', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qi_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QI', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qr_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QR', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qg_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QG', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   qs_crm = load_reshape(spcamCaser, it-ilag, 'CRM_QS', [1 1 1], [inf inf inf], nlon, nlat, nx, nz);
   spcamtstep = strrep(strrep(spcamCaser.name(it-ilag),'spcam_actual_m2005_f09f09_branch.cam.r.',''),'.nc',''); 
   prec_crm = ncread([spcamCaseh1.dir 'spcam_actual_m2005_f09f09_branch_4CRMPREC.cam.h1.' spcamtstep{1} '.nc'],...
                     'CRM_PREC', [1 1 1 1 1], [inf inf inf inf inf ])*3.6e6; %  1 m/s = 100*3600 mm/hr
                      % replace string with blank and obtain the timestep string
   %cape_cam = ncread([spcamCaser.dir 'cape-cin/spcam_actual_m2005_f09f09_branch.cam.rh0.' spcamtstep{1} '-cape-cin.nc'], 'CAPE', [1 1], [inf inf]); 
   %cin_cam  = ncread([spcamCaser.dir 'cape-cin/spcam_actual_m2005_f09f09_branch.cam.rh0.' spcamtstep{1} '-cape-cin.nc'], 'CIN',  [1 1], [inf inf]);
   %cape_crm = ncread([spcamCaser.dir 'cape-cin_crmave/spcam.m2005.crmave.r.' spcamtstep{1} '-cape-cin.nc'],              'CAPE', [1 1], [inf inf]);                  
   %cin_crm  = ncread([spcamCaser.dir 'cape-cin_crmave/spcam.m2005.crmave.r.' spcamtstep{1} '-cape-cin.nc'],              'CIN',  [1 1], [inf inf]);                  
   %%%%%%%%%%%%%%%%%%%%%%%
   % Load LS Variables  %
   %%%%%%%%%%%%%%%%%%%%%%%
   % CAM: [lat lon lev time] CAMO: [lon lat lev time]
   % loadvar(VNAME, CASEINFO, varScalingFactor, numMinusTime, startIndices, numIndices)
   z3_cam     = flipdim(loadvar('Z3',     camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % flip the level dim
                             % geopotential height \int g dz / g_0, flip the height to increase from surface to toa
   %lcwat_cam  = flipdim(loadvar('LCWAT',    camCaser,   1,     0, [1 1 1 it], [inf inf inf 1]), 3); % cloud liquid+ice mixing ratio [kg/kg], 
                                                                                                     % flip the pressure dim for interpolation
   u_cam      = flipdim(loadvar('U',      camCaser,     1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % vel in CAM 
   T_cam      = flipdim(loadvar('T',      camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % temp in CAM [K]
   Q_cam      = flipdim(loadvar('Q',      camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % specific humidity in CAM [kg/kg]
   precc_cam  =         loadvar('PRECC',  camCaserh0,   3.6e6,     0, [1 1 it-ilag],   [inf inf 1]); % convective prec rate [m/s]->[mm/hr]
   precl_cam  =         loadvar('PRECL',  camCaserh0,   3.6e6,     0, [1 1 it-ilag],   [inf inf 1]); % stratiform prec rate [m/s]->[mm/hr]
   prect_cam  =         loadvar('PRECT',  camCaserh0,   3.6e6,     0, [1 1 it-ilag],   [inf inf 1]); % convective prec rate [m/s]->[mm/hr]
   prec_spcam =         loadvar('SPPFLX', spcamCaserh0, 3.6e6,     0, [1 1 30 it-ilag],   [inf inf 1 1]); % sp prec rate [m/s]->[mm/hr]
   dtcond_cam = flipdim(loadvar('DTCOND', camCaserh0,   86400, 0, [1 1 1 it-ilag], [inf inf inf 1]),3); % heating in CAM [K/day]
   spdt       = flipdim(loadvar('SPDT',   spcamCaserh0, 86400, 0, [1 1 1 it-ilag], [inf inf inf 1]),3); % heating in SPCAM [K/day]
%   relhum_tmp = flipdim(loadvar('RELHUM', camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]), 3);    % [%]
   icimr_cam  = flipdim(loadvar('ICIMR',  camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % cloud ice mixing ratio [kg/kg], 
   icwmr_cam  = flipdim(loadvar('ICWMR',  camCaserh0,   1,     0, [1 1 1 it-ilag], [inf inf inf 1]),3); % cloud liquid mixing ratio [kg/kg], 
   lcwat_cam  = icimr_cam + icwmr_cam; clear icimr_cam icwmr_cam;
   display(['finish loading all variables for ilag = ' num2str(ilag)])
   if (itt == 1) % first loop 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  Select Large PC1 lon-lat index  %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [iillt_zm_PC, iCAMll] = getLargePClonlat(it, iillt_zm, 1, bsign, pc_spz, threshold_bigPC); %iCAMll: col-index of iillt_zm on big PC 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Find MCS centroid lon-lat %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [mcs, C] = getMCSlonlat(T_crm, qT_crm, qc_crm, lon, lat, landfrac, noLand, iillt_zm_PC, iCAMll, nlon, nlons, nlons_half, p0, ps, hyai, hybi, hyam, hybm, nz);
      for ic = 1:mcs.ncentroids
         ilons(:,ic) = getlongs(mcs.illcentroids(ic,1),nlon,nlons_half,mcs.illcentroids(ic,2),landfrac, noLand);
      end
      %%%%%%%%%%%%%%%%%%%%%%%
      % MCS heating modes   %
      %%%%%%%%%%%%%%%%%%%%%%% 
      for illct = 1:size(mcs.ilonlat,1)
         spdt_mcs(illct,:)   = squeeze(spdt(mcs.ilonlat(illct,1), mcs.ilonlat(illct,2),:));
         dtcond_mcs(illct,:) = squeeze(dtcond_cam(mcs.ilonlat(illct,1), mcs.ilonlat(illct,2),:));
      end
      [svec_mcs_sp sval_mcs_sp pc_mcs_sp] = svdVar(spdt_mcs'); 
%      coeff = pca(spdt_mcs);
%      coeff 
%      svec_mcs_sp
%      pc_mcs_sp
%pause
      [svec_mcs    sval_mcs    pc_mcs   ] = svdVar(dtcond_mcs');
      %scatter(mcs.llcentroids(:,1),mcs.llcentroids(:,2)); pause
   end % if(itt==1) first loop
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Get neighboring fields of MCS centroid  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for illct = 1:size(mcs.illcentroids,1); 
      ilatc = mcs.illcentroids(illct,2); % MCS center ilat
      %ilons(:,illct) = getlongs(mcs.illcentroids(illct,1),nlon,nlons_half,ilatc,landfrac,noLand);
      for ilonct = 1:nlons
         ilon = ilons(ilonct,illct); % MCS surrounding ilons at fixed ilatc
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %      interpolate variable (p to z)     %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         [pi, pm, dpi, dpm]  = hybrid2p(p0, ps(ilon,ilatc), hyai, hybi, hyam, hybm, nz);
         z0 = phis(ilon,ilatc)/g; % surface geopotential height
         for ix = 1:nx
            for iz = 1:nz
               [rho(iz), ~, rh_crm(ix,iz), qv(ix,iz)]  = density_temp(T_crm(ilon,ilatc,ix,iz), pm(iz),...
                              qT_crm(ilon,ilatc,ix,iz), qc_crm(ilon,ilatc,ix,iz), [], compset); % midpoint density
            end
            [z, dzi] = p2z(rho,dpi,dpm,g,z0,nz);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Interpolate CRM variable p2m %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            varCRM2D(ilonct,ix,1:nzi,itlag,illct,:) = interp1qr(z', ...
                                                         [squeeze(u_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(w_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(T_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(qc_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(qi_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(qr_crm(ilon,ilatc,ix,1:nz)), ...
                                                          squeeze(rh_crm(ix,1:nz))', ...
                                                          ], zint);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             buoyancy freq                   %
            % (bring a parcel from lower to higher level) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            n2 = - g^2 ./ rho(1:nz-1) .* (rho(2:nz)-rho(1:nz-1))./(z(2:nz)-z(1:nz-1)); 
                          % assume n2 misses the highest level, so it'll be at z(1) to z(nz-1)
%            n2crm(ix,1:nzi,itlag,illct) = interp1qr(z(1:nz-1)', n2', zint);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reflectivity dBz Calc          %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            [ WRF_DBZ(ix,:), WRF_DBZ_MAX(ix), Ze(ix,:) ] = wrf_dbz(pm',squeeze(T_crm(ilon,ilatc,ix,:)),...
%                                                        squeeze(qr_crm(ilon,ilatc,ix,:)),...
%                                                        squeeze(qg_crm(ilon,ilatc,ix,:)),...
%                                                        squeeze(qs_crm(ilon,ilatc,ix,:)),...
%                                                        qv_crm(ix,:)');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Convective-Stratiform Separation (inaccurate!!!) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %conv_strat_separation(z,squeeze(w(ilon,ilatc,ix,:)),squeeze(qc(ilon,ilatc,ix,:)),WRF_DBZ_MAX(ix));  
            %[cloudtype, rainrate] = conv_strat_separation();
            [varCRM1D(ilonct,ix,itlag,illct,1), varCRM1D(ilonct,ix,itlag,illct,2)] = conv_strat_separation( ...
                                                         squeeze(prec_crm(ilon,ilatc,ix)), ...
                                                         squeeze(w_crm(ilon,ilatc,ix,1:nz)), ...
                                                         squeeze(qc_crm(ilon,ilatc,ix,1:nz)), ...
                                                         [],[],[],[],[],[],[]); 
         end % ix
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % LS Conv and Strat Rainrate & Frequency    %
         % rainrate_conv(1) rainrate_strat(2)        % 
         % convrain_freq(3) stratrain_freq(4)        % 
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         varCRM0D(ilonct,itlag,illct,1:4) = ...
            [ nansum(varCRM1D(ilonct,(varCRM1D(ilonct,1:nx,itlag,illct,1) == 1),itlag,illct,2))/nx,... 
              nansum(varCRM1D(ilonct,(varCRM1D(ilonct,1:nx,itlag,illct,1) == 0),itlag,illct,2))/nx,...
              nansum(varCRM1D(ilonct,:,itlag,illct,1) == 1)/nx,...
              nansum(varCRM1D(ilonct,:,itlag,illct,1) == 0)/nx ]; 
         varCRM0D(ilonct,itlag,illct,5) = prec_spcam(ilon,ilatc);
         varCAM0D(ilonct,itlag,illct,1:2) = [precc_cam(ilon,ilatc), precl_cam(ilon,ilatc)]; % conv prec, strat prec
         varCAM0D(ilonct,itlag,illct,3) = precc_cam(ilon,ilatc) + precl_cam(ilon,ilatc); % prect = precc+precl
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Stability Index: CAPE, CIN, LI %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         stabilityIndex(ilonct,itlag,illct,1:4) = [cape_cam(ilon,ilatc),cape_crm(ilon,ilatc),...
%                                                     cin_cam(ilon,ilatc),cin_crm(ilon,ilatc)];
%         stabilityIndex(ilonct,itlag,illct,5) = liftedIndex(T_cam(ilon,ilatc,iLI), pm(iLI), ...
%                              1/(1/Q_cam(ilon,ilatc,iLI)-1), (pm(2:end)), (squeeze(T_cam(ilon,ilatc,1:nz)))); 
%                                                      % CAM; use p(2) since LI is calc using surface + 50mb
         % use the CRM bottom 32 points as 32 parcels and
         % lifted to T_env of CAM and CRM
         Tpar  = squeeze(T_crm(ilon,ilatc,:,iLI)); 
         qvpar = squeeze(qT_crm(ilon,ilatc,:,iLI)-qc_crm(ilon,ilatc,:,iLI));
         [stabilityIndex(ilonct,itlag,illct,2), tp500] = liftedIndex( Tpar, pm(iLI), qvpar, ... 
                          pm(2:end), squeeze(mean(T_crm(ilon,ilatc,:,1:nz),3)) ); % CRM, tp500 is the most unstable parcel 
         stabilityIndex(ilonct,itlag,illct,1) = liftedIndex( [], [], [], ...
                          pm(2:end), squeeze(T_cam(ilon,ilatc,1:nz)), tp500 ); % CAM, use the most unstable tp500 from CRM  
%         dps=abs(ps(ilon,ilatc)-5000-p(1:5)); % p-diff
%         iLI = find(dps == min(dps)); % find the pres index of the parcel at 5000Pa above surface pressure 
%         stabilityIndex(ilonct,itlag,illct,5) = liftedIndex(T_cam(ilon,ilatc,iLI), p(iLI), ...
%                         1/(1/Q_cam(ilon,ilatc,iLI)-1), (p(2:end)), (squeeze(T_cam(ilon,ilatc,1:nz)))); 
%                                                      % CAM; use p(2) since LI is calc using surface + 50mb
%         stabilityIndex(ilonct,itlag,illct,6) = liftedIndex(mean(T(ilon,ilatc,:,iLI)), p(iLI), ...
%                         mean(qT(ilon,ilatc,:,iLI)-qc(ilon,ilatc,:,iLI)), (p(2:end)), (squeeze(mean(T(ilon,ilatc,:,1:nz),3)))); % CRM
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Interpolate CAM fields p2z %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         for iz=1:nz
            relhum_cam(iz) = specific2relhum(T_cam(ilon,ilatc,iz),pm(iz),Q_cam(ilon,ilatc,iz));
         end
         varCAM1D(ilonct,1:nzi,itlag,illct,:) = interp1qr(squeeze(z3_cam(ilon,ilatc,1:nz)), ...
                                                                [squeeze(lcwat_cam(ilon,ilatc,1:nz)), ...
                                                                 squeeze(dtcond_cam(ilon,ilatc,1:nz)), ...
                                                                 squeeze(spdt(ilon,ilatc,1:nz)), ...
                                                                 relhum_cam', ...
                                                                 squeeze(u_cam(ilon,ilatc,1:nz)), ...
                                                                 ], zint);
%                                                                 squeeze(relhum_tmp(ilon,ilatc,:)), ...
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % wind shear between  zint(1) and z3km %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         stabilityIndex(ilonct,itlag,illct,3) = abs(varCAM1D(ilonct,i3km,itlag,illct,5) - varCAM1D(ilonct,1,itlag,illct,5))/dz3km ; % |du|/dz
         stabilityIndex(ilonct,itlag,illct,4) = abs(mean(varCRM2D(ilonct,:,i3km,itlag,illct,1),2) - mean(varCRM2D(ilonct,:,1,itlag,illct,1),2))/dz3km ;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  specific mass divergence  %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if (dodivvort)
            [ div(1:nx,1:nzi,itlag,illct), vort(1:nx,1:nzi,itlag,illct) ] = div_vort(varCRM2D(:,:,itlag,illct,1),varCRM2D(:,:,itlag,illct,2),dx,dzi);
         end         
         %%%%%%%%%%%%%%%%%
         %  ke spectrum  %
         %%%%%%%%%%%%%%%%%
         if (dofft)
            FFTke(1:nx/2+1,1:nzi/2+1,itlag,illct) = fftke(varCRM2D(:,:,itlag,illct,1),varCRM2D(:,:,itlag,illct,2),nx,nzi);
         end
      end % ilon
   end % ilatc
   %contourf(reshape(permute(squeeze(qrcam(:,:,itlag,illct,:)),[1,3,2]),nx*(ilonct-1),nz)','linestyle','none'); pause
toc
end % ilag
%
%%%%%%%%%%%%%%%%%%%%
% Combine all parm %
%%%%%%%%%%%%%%%%%%%%
%reflectivity.WRF_DBZ=WRF_DBZ;
%reflectivity.WRF_DBZ_MAX=WRF_DBZ_MAX;
%parm.reflectivity = {'WRF_DBZ','WRF_DBZ_MAX','Ze'}
%reflectivity.Ze=Ze;
parm.varCRM2D = {'u','w','T','qc','qi','qr','rh'}; % qli = qc+qi = qT-qv, these are variables for CRM
parm.varCRM1D = {'cloudtype','rainrate'}; % cloudtype: convective=1 stratiform=0, rainrate [mm/hr]
parm.varCRM0D = {'precc','precl','freq_precc','freq_precl','prect'}; % cloudtype: convective=1 stratiform=0, rainrate [mm/hr]
parm.varCAM1D = {'lcwat','dtcond','spdt','relhum','u'}; % lcwat = qc+qi 
parm.varCAM0D = {'precc','precl','prect'}; 
%parm.stabilityIndex= {'cape_cam','cape_crm','cin_cam','cin_crm','LI_cam','LI_crm','du3km_cam','du3km_crm'};
parm.stabilityIndex= {'LI_cam','LI_crm','du3km_cam','du3km_crm'};
parm.lon = lon; parm.lat = lat;
parm.hyam = hyam; parm.hybm = hybm; parm.hyai = hyai; parm.hybi = hybi;
parm.phis = phis; parm.p0 = p0; 
parm.nx = nx; parm.zint = zint; parm.dzi=dzi;
parm.nzi=nzi; parm.nlag=nlag; parm.nlons=nlons;
parm.lagt=lagt; parm.nlagt = nlagt; parm.dx=dx;
parm.iillt_zm_PC=iillt_zm_PC; 
parm.mcs = mcs;
%parm.iCAMll_ocn = iCAMll_ocn;



%%%%%%%%%%%%%%%%%%%%%%%
%   Save variables    %
%%%%%%%%%%%%%%%%%%%%%%%
save(fname,'parm','varCAM0D','varCAM1D','varCRM0D','varCRM1D','varCRM2D','stabilityIndex','mcs','svec_mcs_sp','sval_mcs_sp','pc_mcs_sp', 'svec_mcs', 'sval_mcs', 'pc_mcs')
%save(fname,'parm','varCAM0D','varCAM1D','varCRM0D','varCRM1D','varCRM2D','stabilityIndex','mcs')
%save(fname,'parm','var','n2crm','div','vort','FFTke','qr')
%save(['T_qc_qt_qr_PC1_' num2str(it) '.mat'],'parm','Tcrm','qccrm','qTcrm','qrcrm')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iillt_zm_PC iCAMll] = getLargePClonlat(it, iillt_zm, iEOF, bsign, pc_spz, threshold_bigPC)
   % iillt_zm:    matrix for lon-lat-time indices associated to ZM deep
   %              convection
   % iillt_zm_PC: matrix of iillt_zm conditioned on big iEOF-PC 
   % iCAMll:      column indices of iillt_zm_PC at it-timestep
   for ib = 1:4 % basin index 
      iitmp = iillt_zm{ib}(:,find(bsign(ib)*pc_spz{ib}(:,iEOF) > threshold_bigPC )); 
                          % pick out the big PC1 amplitude lon, lat, time indices 
                          % and use it to study the CRM KE spectrum
                          % rows: 1(lon) 2(lat) 3(time)
      iillt_zm_PC{ib} = sortrows(iitmp',3)'; % sort the rows according to the 
                          % 3rd (time) column of iitmp'<--- this sorting might not be necessary
      iCAMll{ib} = find(iillt_zm_PC{ib}(3,:)==it); % find the column indices of 
                          % iillt_zm_PC{ib} for lat & lon at it-timestep, 
                          % use these lat & lon for all lag-lead timesteps
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


function [mcs, C] = getMCSlonlat(T,qT,qc,lon,lat,landfrac,noLand,iillt_zm_PC,iCAMll,nlon, nlons, nlons_half, p0, ps, hyai, hybi, hyam, hybm, nz)
   % Get MCS cluster center lon-lat and ilon-ilat
   cond_plot=0;  % scatter plot for the MCS points and cluster center points
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Calc LI for all big iEOF-PC points %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for ib = 1:4
      if (isempty(iCAMll{ib})==0) % if not empty 
         illct = 1; % counter for the lon lat indices with zonal lons over ocn only
         for ill = iCAMll{ib}; % lon lat column index of iillt_zm_PC{ib}
            ilon = iillt_zm_PC{ib}(1,ill);
            ilat = iillt_zm_PC{ib}(2,ill);
            ilons = getlongs(ilon,nlon,nlons_half,ilat,landfrac,noLand); % obtain the zonal 
                          % longitudes surrounding the large PC lon, at a fixed latitude, without land 
            if (ilons == 0); % if there are land covered, skip/discard this lon-lat index
               continue
            end
            iCAMll_ocn{ib}(illct) = ill; % column indices of iillt_zm_PC{ib} 
                       % that has big PC1 lon extending zonally over ocean only.
            for iln=1:nlons 
               ilon = ilons(iln);
               iLI=1;
               [~, pm, ~, ~]  = hybrid2p(p0, ps(ilon,ilat),hyai,hybi,hyam,hybm,nz);
               Tpar  = squeeze(T(ilon,ilat,:,iLI)); 
               qvpar = squeeze(qT(ilon,ilat,:,iLI)-qc(ilon,ilat,:,iLI));
               [LI{ib}(iln,illct), ~] = liftedIndex( Tpar, pm(iLI), qvpar, ... 
                                   pm(2:end), squeeze(mean(T(ilon,ilat,:,1:nz),3)) ); 
                       % CRM, tp500 is the most unstable parcel 
            end % iln
            illct = illct+1;
         end % ill 
      end % if(isempty(iCAMll{ib})==0)
   end % ib
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Find MCS ilon-ilat        %
   % (big PC points with LI<0) % 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % get mcs.lonlat mcs.ilonlat, mcs.centroids
   ilon=[];ilat=[];
   for ib=1:4
      iilltocn=[];
      for iilon=1:nlons
         ii1 = find(squeeze(LI{ib}(iilon,:))<0); % LI is from MCS index
         iilltocn = [ii1,iilltocn]; % indices of iCAMll_ocn{ib} 
      end
      if (isempty(iilltocn))
         continue 
      end 
      ill  = iCAMll_ocn{ib}(iilltocn); % column indices of iillt_zm_PC{ib}
      ilon = cat(2,iillt_zm_PC{ib}(1,ill),ilon); % concat all basins to one ocean
      ilat = cat(2,iillt_zm_PC{ib}(2,ill),ilat);
   end
   mcs.lonlat  = unique([lon(ilon),lat(ilat)],'rows'); % remove all the overlapped lon-lat pairs
   mcs.ilonlat = unique([ilon',ilat'],'rows');
   cutoff=2.8; % used for 'euclidean' 'single' distance
   [C Z] = myclusterdata(mcs.lonlat,'euclidean','single','distance',cutoff); % C: cluster index array; Z: pairwise-linkage-distance matrix
   ii=1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Find MCS cluster central ilon-ilat %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for ic = 1:max(C)% index of clusters 
      idc = find(C==ic); % find the mcs.lonlat row indices corresponding to each cluster 
      if (numel(idc)<3); continue; end % skip clusters with less than 3 grid points
      [~, mcs.llcentroids(ii,1:2)] = kmedoids(mcs.lonlat(idc,:),1,'algorithm','small','distance','euclidean'); % obtain all the center points in each MCS system 
      mcs.illcentroids(ii,1) = find(lon==mcs.llcentroids(ii,1));
      mcs.illcentroids(ii,2) = find(lat==mcs.llcentroids(ii,2));
      mcs.units(ii) = numel(idc); % # of units in a cluster
%      mcs.llcentroids(ii,1:2) = median(mcs.lonlat(idc,:)); % obtain all the center points in ech MCS system 
%      mcslon = (abs(lon-mcs.llcentroids(ii,1)));
%      mcslat = (abs(lat-mcs.llcentroids(ii,2))); % if llcentroid=0, then it'll give two lat values
%      mcs.illcentroids(ii,1) = min(find(mcslon==min(mcslon))); 
%      mcs.illcentroids(ii,2) = min(find(mcslat==min(mcslat)));% might be more than two 
      ii=ii+1; % counter for greater than 3 grid point clusters (defined as MCS systems)
   end
   mcs.ncentroids=size(mcs.illcentroids,1); % total # of centroids
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % test plot how the cluster looks like on a lon-lat map %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (cond_plot)
      figure;
      for i=1:size(C,1);
         scatter(mcs.lonlat(i,1),mcs.lonlat(i,2),'b'); hold on; 
         h(i)=text(mcs.lonlat(i,1),mcs.lonlat(i,2),num2str(C(i))); set(h(i),'FontSize',10) % label the clusters in number
      end; %pause
      scatter(mcs.llcentroids(:,1),mcs.llcentroids(:,2),'r'); hold on; pause
      % use the T cluster index as MCS cluster index
   end

function [var]       = load_reshape(Case, it, varname,  start, count, nlon, nlat, nx, nz); 
   var  = ncread([Case.dir Case.name{it}],varname,start,count); 
   var  = reshape(var,nlon,nlat,nx,nz);

function [div, vort] = div_vort(u, w, dx, dzi)
   [dudx dudz] = gradient(u,dx,dzi);
   [dwdx dwdz] = gradient(w,dx,dzi);
   div  = dudx + dwdz;
   vort = dudz - dwdx;

function [FFTke]     = fftke(u, w, nx, nzi)
   FFTu = fft2(u); % fft: columnwise zonal fft for each pressure level independently
   FFTw = fft2(w); % fft2: row and column-wise fft, i.e., fft(fft(uw))
   FFTmagu = abs([FFTu(1,:); 2*FFTu(2:nx/2,:); FFTu(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagw = abs([FFTw(1,:); 2*FFTw(2:nx/2,:); FFTw(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
   FFTmagu = ([FFTmagu(:,1), 2*FFTmagu(:,2:nzi/2), FFTmagu(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTmagw = ([FFTmagw(:,1), 2*FFTmagw(:,2:nzi/2), FFTmagw(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
   FFTke = 0.5*(FFTmagu.^2 + FFTmagw.^2);

function ilons     = getlongs(ilon, nlon, nlons_half, ilat, landfrac, noLand)
   ilons = ilon-nlons_half:ilon+nlons_half;
   if (ilons(1)<=0);        ilons(find(ilons<=0))      = [nlon-sum(ilons<=0)+1:nlon]; end 
                 % if long fell out the left/right boundary, use the right/left boundary periodically
   if (ilons(end)>=nlon+1); ilons(find(ilons>=nlon+1)) = [1:sum(ilons>=nlon+1)];      end
   if (noLand) % if true, discard land points if landfrac is greater than zero
      if (sum(landfrac(ilons,ilat))~=0) % if there are land indices involved, discard this sample
         ilons = 0;
      end
   end

function [RH] = specific2relhum(T,P,q)
   qv = q/(1-q);
   qvs = r_sub_s(P,T);
   RH = qv/qvs*100;


function [rho T_rho RH qv] = density_temp(T,P,qT,qc,qi,compset)
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
   if (compset == 'm2005')
      qv = qT - qc;
      %   qvs = r_sub_s(P,T); % saturation mixing ratio, used for
      %   calculating RH = qv/qvs
      %end
   elseif (compset == 'sam1mom')
      qv  = qT - qc - qi; % qt = total water; qc = cloud liquid; qi = cloud ice; qv = cloud vapor
   end
   if (nargout >= 3) % even if T_rho is ~ not called it will be counted as nargout
      qvs = r_sub_s(P,T); % saturation mixing ratio, used for calculating RH = qv/qvs, r_sub_s is used in LI function
      RH  = qv/qvs*100;
   end
   T_rho = T*(1+zvir*qv);
   %T_rho = T*(1+qv/0.622)/(1+qv+qc+qi);
   %T_rho1 = T.*(1+0.61*qv - qc - qi); % approximate density temperature
   rho = P.*rRd./T_rho ; % P = rho Ra T_rho => rho = P * 1/Ra * 1/T_rho (1/Rd=0.003484)

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
   dp     = abs(p-50000); % find the closest p-lev to 50000
   k500   = find(dp == min(dp));
   te500  = t(k500) + (t(k500+1)-t(k500))/(p(k500+1)-p(k500)) * (50000-p(k500)); % env temp interpolated to closest k500 level
   if (nargin<6) % compare all parcels and find the most unstable (max(thelcl)) one and output it
      parm1  = (P0/pp)^RCP;
      for ip = 1:length(Tp) % number of parcels
         th     = Tp(ip) * parm1; % potential temp of the parcel 
         Td     = tdew(pp,qv(ip)); % dew point temp of the parcel
         tlcl   = t_lcl(Tp(ip),Td); % temp at lcl 
         plcl   = P0 * (tlcl/th)^CPR; % pressure at lcl
         thelcl(ip) = theta_e(plcl,tlcl,qv(ip),tlcl);
      end
      imaxp = find(thelcl == max(thelcl)); % max thelcl parcel index
      thelcl = thelcl(imaxp(1)); % find the most unstable parcel, if more than one max, just use the first one since all values are the same
      tp500  = compT_fr_The(thelcl,50000); % temp of the parcel at p500
   end
   LI  = te500 - tp500; % env - parcel 

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


function [cloudtype rainrate] = conv_strat_separation(prec, w, qc, qr, dbz, dbz_max, Ze, T, z, dzi)
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
      elseif ( prec > 0.1 ) %[mm/hr]
         cloudtype = 0;
      else
         cloudtype = NaN;
      end
   end
      
   % stratiform:
   % if local point precip
   % > 0.1 mm/hr 
   % rainfall rate and reflectivity uses 
   % Z = a * R^b 
   % where Z [mm^6 m^-3] is the reflectivity, and R [mm h^-1] is the
   % rainfall rate, for Marshall-Palmer relation a=200, b=1.5
   % qc_thres = 0.0005; % [kg/kg]
   % w_thres = 3; %[m/s]
   % if ( qc > qc_thres) | w > w_thres )
   % RdRv = 0.622; 
   %   qv = qv+1e-8; 
   %   e = p*qv/(qv+RdRv);
   %   pd = p - e;
   %   rhoair = pd/(Rd*T);
   %   rhoh2o = qr*rhoair; %[kg/m3]

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % convective-stratiform separation criteria by radar reflectivity  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % http://www.wdtb.noaa.gov/courses/MRMS/ProductGuide/PrecipQualityControls/conv-strat-precip-sep.php
   % 1. max dBZ (<25km) >35 dBZ (since >45 is the higher threshold for convective, but
   %    actually only need 35)
   % 2. rainrate calc use Marshall-Palmer formula a=200, b=1.5
   %    https://en.wikipedia.org/wiki/DBZ_(meteorology)
   %    Z = a * R^b, dBZ \propto 10*log_10(Z/Z0)  
%   TC = T - 273.15; % temp in celcius
%   iz13 = find(abs(z-1300)  == min(abs(z-1300)));  % 1.3km height index
%   iz25 = find(abs(z-25000) == min(abs(z-25000)));
%   iz0C = find(abs(TC)      == min(abs(TC))); % 0C freezing level
%   % below 25km
%   if (any(qr(1:iz0C) > 0.01*1e-3)) % if rain mixing ratio > 0.01g/kg=0.01*1e-3kg/kg below freezing level
%      if (z(iz0C) < 2000) % freezing level < 2km
%         cloudtype=0; % stratiform
%      else % freezing level >= 2km
%         if (max(dbz(1:iz13)) > 45 | max(w(1:iz13))>3 | max(qc(1:iz13))>0.5*1e-3 ) 
%                                         % if max reflectivity below 1.3km > 45dBZ
%                                         % if max w below 1.3km > 3m/s 
%                                         % if max qc below 1.3km > 0.5g/kg=0.5*1e-3kg/kg 
%            cloudtype=1; % convective
%         else
%            cloudtype=0;
%         end
%      end
%   else % above 25km and below 250km (unfortunately spcam doesn't reach above 28km, so can't use this part)
%      if ( any(dbz(iz25:end)>0) )
%         vilwc = vil(Ze(iz25:end),dbz(iz25:end),dzi(iz25:end)); % vertically integrated liquid water content [kg/m^2]
%         vilwc = vil(Ze,dbz,dzi); % vertically integrated liquid water content [kg/m^2]
%      else
%         cloudtype=NaN;
%      end
%      if (exist('cloudtype')==0) % if cloudtype nonexist 
%         if (vilwc < 6.5)
%            cloudtype=0; % stratiform
%         else % vilwc >=6.5
%            izN10C = find(abs(TC+10) ==  min(abs(TC+10))) % -10C level
%            if ( dbz(izN10C) < 30 ) % if -10C level reflectivity < 30 dBZ
%               cloudtype=0; % stratiform
%            else % if -10C level reflectivity >= 30 dBZ 
%               cloudtype=1; % convective
%            end   
%         end
%      end
%   else
%      cloudtype=NaN;
%   end
%   if (dbz_max == -30)
%      rainrate = 0;
%   else
%      rainrate = (10^(dbz_max/10)/200)^(5/8); % [mm/hr] 8/5=1.6
%   end



   
function vilwc = vil(Ze,dbz,dzi)
   % Greene and Clark 1972
   % http://www.wdtb.noaa.gov/courses/MRMS/ProductGuide/SevereWeather/vertically-integrated-liquid.php
   % VIL:=vertically integrated liquid-water content [kg/m^2] derived
   % from radar reflectivity ( check consistency with actual LWC
   % Ze [mm^6 m^-3]
   % dbz = 10 * log10(Z/Z0)
   % VIL [kg m^-2]
%     dbz(dbz>56)=56; % to exclude contribuitions from ice, set to 56 dBZ
     ipos = find(dbz>0);
     Ze(dbz>56)=10^(56/10); % to exclude contribuitions from ice, set to 56 dBZ
     vilwc = 3.44*1e-6*dot(Ze(ipos).^(4/7),dzi(ipos));
%     vilwc = 3.44*1e-6*dot(dbz(ipos).^(4/7),dzi(ipos));
