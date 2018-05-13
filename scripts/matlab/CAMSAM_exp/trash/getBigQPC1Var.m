function var = getBigQPC1Var(varname)
% Purpose: obtain the variables on the space-time index exceeding Q 
% big PC1 threshold (> 5K/day).

run /nethome/gchen/scripts/matlab/startup.m
Case.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
Caserh0.dir    = Case.dir;
%Caserh0.dir = '/projects/rsmas/kirtman/gchen/archive/spcam_cam_rh0_m2005/atm/rest';
dir_r       = dir([Case.dir '/spcam_actual_m2005_f09f09_branch.cam.r.*.nc']);
dir_rh0     = dir([Caserh0.dir '/spcam_actual_m2005_f09f09_branch.cam.rh0.*.nc']);
Case.name   = {dir_r.name};
Caserh0.name   = {dir_rh0.name};
CaseDim.name = 'spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-01800' 
CaseDim.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
lat  = loaddim('lat',CaseDim);
lon  = loaddim('lon',CaseDim);
lev  = loaddim('lev',CaseDim);  % level centers
ilev = loaddim('ilev',CaseDim); % level interfaces, two interfaces surround one level center
                                % index 1 (TOA pressure) to nlev (surface pressure) 
dlev = ilev(2:end)-ilev(1:end-1); % level thickness (>0)
nlat=numel(lat);
nlon=numel(lon);
nlev=numel(lev);
nx    = 32; % evenly spaced
km2m  = 1000; % 1km = 1000m
dx    = 4*km2m; % 4 km grid box
nz    = 28; % levels of CAM from surface to nz levels (according to cam/crm_physics.F90 line 992)
nlag  = 7; % extract the nlag timesteps for each grid point that exceeds the PC1 threshold at 0-lag timestep
lagind = nlag:-1:-nlag;
% dtlag = 30; % 1800 sec = 30 min for CAM timestep
% timelag = nlag*dtlag:-dtlag:0; % 0 lag is the big PC1 event timestep
g = 9.8; 
crmlev  = lev(nlev:-1:nlev-nz+1); % 1 (surface pressure) to nlev (TOA)
crmdp(1) = crmlev(1)-ilev(end); % p-diff between 1st level to bottom-most pressure
crmdp(2:nz)= crmlev(2:nz) - crmlev(1:nz-1); % the center to center pressure thickness of each level
crmdlev = dlev(nlev:-1:nlev-nz+1);
%
%zint = linspace(1000,31000,16)'; % interpolated height with 2km grid-spacing
nzi  = 180; % # of interpolated z levels, nzi=61(~500m), nzi=180(~150m)
maxz = 28000; % ~31500m are most of the heighest heights of each vert columns
minz = 100; % ~ 55m are most of the lowest heights of each vert columns
zint = linspace(minz,maxz,nzi)'; % interpolated height with ~500m grid-spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Heating profile SVD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
threshold_bigPC = 10; % the 1st EOF mode largest magnitude vertically is 0.5, and the average PC has a vertical maximum of 5 K/day heating rate, so a threshold of 10 to get 10*0.5=5 is necessary
iEOF = 1; % EOF1 index
nzz=17;
start = [1 1 1]; % starting index for lon(1:288) lat(1:192) nx*nz(1:896) (maybe use lat(76:117) 20 degrees north-south?) 
count = [inf inf inf]; % number of index to read in var, inf means use all elements in that dimension
%      FFTke{ib} = zeros(nx/2+1,nzi/2+1,length(lagind),length(iCAMll));
for ilag = lagind % lagging big PC1's index by ilag, so ilag=0 is the big PC1 index
   itlag = nlag + 1 - ilag; % save nlag at first index, no lag at last index
%         u  = ncread([Case.dir Case.name{it-ilag}],'CRM_U',start,count);  u  = reshape(u,nlon,nlat,nx,nz);% 2D without meridional velocity
%         w  = ncread([Case.dir Case.name{it-ilag}],'CRM_W',start,count);  w  = reshape(w,nlon,nlat,nx,nz); 
         %w(:,nlag-ilag+1) = ncread([Case.dir cases{iillt_zm_PC1{ib}(3,ii)-ilag}],'CRM_W',start,count);
%         T  = ncread([Case.dir Case.name{it-ilag}],'CRM_T',start,count);  T  = reshape(T,nlon,nlat,nx,nz);
%         qc = ncread([Case.dir Case.name{it-ilag}],'CRM_QC',start,count); qc = reshape(qc,nlon,nlat,nx,nz);
%         qT = ncread([Case.dir Case.name{it-ilag}],'CRM_QT',start,count); qT = reshape(qT,nlon,nlat,nx,nz);
%         sppflx = ncread([Caserh0.dir Caserh0.name{it-ilag}],'SPPFLX',start,count); % SPPFLX(lev,lat,lon) total of 30 lev % p-sfc starts from last index
   if (ilag == nlag) % nlag is the first loop for ilag
      for ib = 1:4 % basin index 
         tic% pick out vertical columns at (ilon, ilat) exceeding EOF1 threshold
         if ib==2 % only WPAC basin (ib=2) has an opposite sign EOF1, the rest here are all positive sign 
            iitmp = iillt_zm{ib}(:,find(pc_spz{ib}(:,iEOF) < -threshold_bigPC )); % pick out the big PC1 amplitude lon, lat, time indices and use it to study the CRM KE spectrum
         else 
            iitmp = iillt_zm{ib}(:,find(pc_spz{ib}(:,iEOF) > threshold_bigPC )); % rows: 1(lon) 2(lat) 3(time)
         end
         iillt_zm_PC{ib} = sortrows(iitmp',3)'; % sort the rows according to the 3rd (time) column of iitmp'
         iCAMll{ib} = find(iillt_zm_PC{ib}(3,:)==it); % find the column indices for lat & lon at it-timestep, use these lat & lon for all lag-lead timesteps
      end
   end
   for ib = 1:4
      if (isempty(iCAMll{ib})==0) % if not empty 
         ill_CAM = 1; % counter for the lon lat indices for CAM
         for ill = iCAMll{ib}; % lon lat index 
            ilon = iillt_zm_PC{ib}(1,ill);
            ilat = iillt_zm_PC{ib}(2,ill);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % vertical interpolate u & w (p to z) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for ix = 1:nx
               for iz = 1:nz
                  rho(iz) = density_temp(T(ilon,ilat,ix,iz),crmlev(iz),qc(ilon,ilat,ix,iz),qT(ilon,ilat,ix,iz));
                  dz(iz)  = -crmdp(iz)/(rho(iz)*g); %dp = rho * g * dz
                  if (iz==1)
                     z(1) = dz(1);
                  else
                     z(iz) = z(iz-1) + dz(iz);
                  end;
               end;
               N2{ib}(ix,1:nz-1,itlag,ill_CAM) =  g^2 * (rho(1:nz-1)-rho(2:nz))./(crmlev(1:nz-1)-crmlev(2:nz))'; % N = - g/rho * drho/dz 
                                                                                                            %   = - g/rho * dp/dz*drho/dp 
                                                                                                            %   = - g/rho * -rho*g * drho/dp =  g^2 * drho/dp                
               % interpolate u w on geopotential heights to linear heights
               uw(ix,1:nzi,1:2) = interp1qr(z',[squeeze(u(ilon,ilat,ix,1:nz)), squeeze(w(ilon,ilat,ix,1:nz))], zint);
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            % vertical wind shear %
            %%%%%%%%%%%%%%%%%%%%%%%
            tmp(1:nx,1)    = u(ilon,ilat,:,1)/crmdlev(1); % u_z = u(1) - 0 % bottom-most velocity equals zero
            for ix=1:nx
               tmp(ix,2:nz) = squeeze(u(ilon,ilat,ix,2:end) - u(ilon,ilat,ix,1:end-1))./crmdlev(2:nz);
            end
            u_p{ib}(1:nz,itlag,ill_CAM) = mean(tmp,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % precip rate (vert column) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spprec{ib}(:,itlag,ill_CAM) = sppflx(ilon,ilat,:); % save sp-cam precipitation rate 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  FFT u, w for specific kinetic energy spectrum %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            FFTu = fft2(uw(1:nx,1:nzi,1)); % fft: columnwise zonal fft for each pressure level independently
            FFTw = fft2(uw(1:nx,1:nzi,2)); % fft2: row and column-wise fft, i.e., fft(fft(uw))
            FFTmagu = abs([FFTu(1,:); 2*FFTu(2:nx/2,:); FFTu(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
            FFTmagw = abs([FFTw(1,:); 2*FFTw(2:nx/2,:); FFTw(nx/2+1,:)])./nx;% only half + 1 is needed for real number fft
            FFTmagu = ([FFTmagu(:,1), 2*FFTmagu(:,2:nzi/2), FFTmagu(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
            FFTmagw = ([FFTmagw(:,1), 2*FFTmagw(:,2:nzi/2), FFTmagw(:,nzi/2+1)])./nzi;% only half + 1 is needed for real number fft
            FFTke{ib}(1:nx/2+1,1:nzi/2+1,itlag,ill_CAM) = 0.5*(FFTmagu.^2 + FFTmagw.^2); 
%            tind(ill_CAM) = iillt_zm_PC{ib}(3,ill)-ilag; % save nlag+1 time samples before and including the PC1 heating profile exceeding the threshold
         %plot(mean(P1(:,:,iit),2)); % vertically averaged zonal KE spectrum
            ill_CAM = ill_CAM+1; % counter for the lon-lat CAM indices at a fixed time
         end
      end
%      itt = itt+1 % number of sorted timesteps in iitmp, 11:10:500 thus has a total of 49 timesteps
%   save('PC1_fft-basinwise_' num2str(it) '.mat','P1','iillt_zm_PC',)  % the fft's of CRM kinetic energy are done on big PC1 vertical profiles at lat-lon-time indices       
   end
toc
end
