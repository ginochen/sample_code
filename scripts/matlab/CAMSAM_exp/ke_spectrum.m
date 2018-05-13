function [FFTke iillt_zm_PC u_p N2] = ke_spectrum(it)
% Purpose : Calculate the Kinetic Energy Spectrum and Vertical Wind
% Shear from SPCAM embedded
% CRM columns nx*nz=32*28=896
%
% 1) Find the lon,lat,time grid-point with large PC1 EOF1 Deep Convection vertical Heating profile
% 2) Use the time index "it=11:10:500", exclude the it=1, 
%    selected from 1 to 755 timesteps of /projects/rsmas/kirtman/gchen/archive/F_2000_4SPCAM_m200501/atm/rest/F_2000_4SPCAM_m200501.cam.rh0.DTCOND.nc
% pc_spz = U*S <----- this is the amplitude of Heating projected on the EOFs
% svec_spz = V
% data = U*S*V' = pc_spz*svec_spz'
% find the indices of big pc values in each basin, and use iillt_zm to
% locate their lat lon time indices

run /nethome/gchen/scripts/matlab/startup.m
% CAUTION!!! spcam_actual_m2005_f09f09_branch.cam.r.*.nc starts from
% 0001-01-14-01800 instead of 0001-01-14-03600 which the rh0 variables
% starts from
Case.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/'
%Case.dir    = '/projects/rsmas/kirtman/gchen/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
dir_r       = dir([Case.dir '/spcam_actual_m2005_f09f09_branch.cam.r.*-*-*-*.nc']);
Case.name   = {dir_r.name};
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
lagind = nlag:-1:-nlag; % actual lag index 
lagind_extra = nlag+1:-1:-nlag; % extra lag is to get the instantaneous field by stracting previous timestep
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
% define vert shear on p-lev
      % define fft variables
%      u   = zeros(nlon,nlat,nx*nz,nlag+1); % don't use zeros(1,nx*nz), its slower to ncread
%      w   = zeros(nlon,nlat,nx*nz,nlag+1); 
      % define p2z interpolation variables
%      T   = zeros(nlon,nlat,nx*nz,nlag+1); 
%      qc  = zeros(nlon,nlat,nx*nz,nlag+1); 
%      qT  = zeros(nlon,nlat,nx*nz,nlag+1); 
%      z   = zeros(nz,1);
%      dz  = zeros(nz,1);
%      rho = zeros(nz,1);
%      uw  = zeros(nx,nz,2); % 2 for each u & w
%      ke = zeros(nx,nzz,nlag+1,length(iCAMll));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Mass for KE 
% large-scale (hydrostatic assumption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m = \rho * V = \rho * dz * dx * dy = \rho * dp/(\rho * g) * dx * dy = dp * dx * dy / g <----- dp is the only varying weight to KE!
%m = crmdlev * dx * dx / g; % the second dx is dy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Heating profile SVD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
threshold_bigPC = 10; % the 1st EOF mode largest magnitude vertically is 0.5, and the average PC has a vertical maximum of 5 K/day heating rate, so a threshold of 10 to get 10*0.5=5 is necessary
iEOF = 1; % EOF1 index
nzz=17;
start = [1 1 1]; % starting index for lon(1:288) lat(1:192) nx*nz(1:896) (maybe use lat(76:117) 20 degrees north-south?) 
count = [inf inf inf]; % number of index to read in var, inf means use all elements in that dimension
%11:10:500 <---- run this in parallel
for ilag = lagind % lagging big PC1's index by ilag, so ilag=0 is the big PC1 index
   itlag = nlag + 1 - ilag; % save nlag at second index, no lag at last index
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %      NCREAD var here     %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   u(:,:,:,itlag)      = ncread([Case.dir Case.name{it-ilag}],'CRM_U',start,count); % u  = reshape(u,nlon,nlat,nx,nz);% 2D without meridional velocity
   w(:,:,:,itlag)      = ncread([Case.dir Case.name{it-ilag}],'CRM_W',start,count); % w  = reshape(w,nlon,nlat,nx,nz); 
  %w(:,nlag-ilag+1) = ncread([Case.dir cases{iillt_zm_PC1{ib}(3,ii)-ilag}],'CRM_W',start,count);
   T(:,:,:,itlag)      = ncread([Case.dir Case.name{it-ilag}],'CRM_T',start,count); % T  = reshape(T,nlon,nlat,nx,nz);
   qc(:,:,:,itlag)     = ncread([Case.dir Case.name{it-ilag}],'CRM_QC',start,count);% qc = reshape(qc,nlon,nlat,nx,nz);
   qT(:,:,:,itlag)     = ncread([Case.dir Case.name{it-ilag}],'CRM_QT',start,count);% qT = reshape(qT,nlon,nlat,nx,nz);
   qr(:,:,:,itlag)     = ncread([Case.dir Case.name{it-ilag}],'CRM_QR',start,count);% rain mixing ratio
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get Current timestep instantaneous field   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u  = u(:,:,:,2:end)  - u(:,:,:,1:end-1); % subtraction gives instantaneous velocity instead of accumulated (for monthly averaging)
%w  = w(:,:,:,2:end)  - w(:,:,:,1:end-1);
%T  = T(:,:,:,2:end)  - T(:,:,:,1:end-1);
%qc = qc(:,:,:,2:end) - qc(:,:,:,1:end-1);
%qT = qT(:,:,:,2:end) - qT(:,:,:,1:end-1);
u  = reshape(u,nlon,nlat,nx,nz,nlag*2+1);
w  = reshape(w,nlon,nlat,nx,nz,nlag*2+1); 
T  = reshape(T,nlon,nlat,nx,nz,nlag*2+1);
qc = reshape(qc,nlon,nlat,nx,nz,nlag*2+1);
qT = reshape(qT,nlon,nlat,nx,nz,nlag*2+1);
%
for ilag=lagind
   itlag = nlag + 1 - ilag;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Large PC1 lon-lat index selection  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %      interpolate variable (p to z)     %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for ix = 1:nx
               for iz = 1:nz
                  rho(iz) = density_temp(T(ilon,ilat,ix,iz,itlag),crmlev(iz),qc(ilon,ilat,ix,iz,itlag),qT(ilon,ilat,ix,iz,itlag));
                  dz(iz)  = -crmdp(iz)/(rho(iz)*g); %dp = rho * g * dz
                  if (iz==1)
                     z(1) = dz(1);
                  else
                     z(iz) = z(iz-1) + dz(iz);
                  end;
               end;
               N2{ib}(ix,1:nz-1,itlag,ill_CAM) =  g^2 * (rho(1:nz-1)-rho(2:nz))./(crmdp(2:nz))'; % N = - g/rho * drho/dz 
                                                                                                            %   = - g/rho * dp/dz*drho/dp 
                                                                                                            %   = - g/rho * -rho*g * drho/dp =  g^2 * drho/dp                
               % interpolate u w on geopotential heights to linear heights
               uw(ix,1:nzi,1:2) = interp1qr(z',[squeeze(u(ilon,ilat,ix,1:nz,itlag)), squeeze(w(ilon,ilat,ix,1:nz,itlag))], zint);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  specific mass divergence  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %for iz=1:nz-1
            %   tmp1(1:nx,iz,ill_CAM)  = (rh0(1:nx,iz+1)*squeeze(w{ib}(ilon,ilat,1:nx,iz+1,itlag)) - rho(1:nx,iz)*squeeze(w{ib}(ilon,ilat,1:nx,iz,itlag)))/crmdp(iz+1);
            %end
            %dRhoW{ib}(1:nz-1,itlag,ill_CAM) = mean(tmp1,1);
            %dRhoU{ib}(1:nz,itlag,ill_CAM) = mean((rh0(2:nx,1:nz).*squeeze(u{ib}(ilon,ilat,2:nx,1:nz,1))-rho(1:nx-1,1:nz).*squeeze(u{ib}(ilon,ilat,1:nx-1,1:nz,1)))/dx,1);
            %%%%%%%%%%%%%%%%%%%%%%%
            % vertical wind shear %
            %%%%%%%%%%%%%%%%%%%%%%%
            %tmp(1:nx,1)    = u(ilon,ilat,:,1,itlag)/crmdlev(1); % u_z = u(1) - 0 % bottom-most velocity equals zero
            %for iz=1:nz-1
            %   tmp2(ix,iz) = squeeze(u(ilon,ilat,1:nx,iz+1,itlag) - u(ilon,ilat,1:nx,iz,itlag))./crmdp(iz+1);
            %end
            %u_p{ib}(1:nz-1,itlag,ill_CAM) = mean(tmp2,1);
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
%         ncid = netcdf.open([Case.dir Case.rname{iillt_zm_PC1{ib}(3,ii)-ilag}]);
%         varid_u = netcdf.inqVarID(ncid,'CRM_U');
%         varid_w = netcdf.inqVarID(ncid,'CRM_W');
%         u = netcdf.getVar(ncid,varid_u,start,count);
%         w = netcdf.getVar(ncid,varid_w,start,count);
%         netcdf.close(ncid)
%contourf(1:32,crmlev(1:nzz),ke(:,:,iit)'); caxis([1e08 5e09])
%set(gca,'YDir','Reverse')
%pause
%         KE = fft(ke); % columnwise zonal fft for each pressure level independently
%         P2 = abs(KE); % get the magnitude for all FFT, the magnitude is the coefficient of the shifted cosine wave from the sum of the sine and cosine wave at the same k, with geometry seen in p162 of Smith99.
%         P1 = P2(1:nx/2+1,:); % only half + 1 is needed for real number fft
%         P1(2:end-1,:) = 2/nx*P1(2:end-1,:);
%         P1(1,:)       = 1/nx*P1(1,:);
%         P1(end,:)     = 1/nx*P1(end,:);
%         tind = iillt_zm_PC1{ib}(3,ii)-ilag; % save nlag+1 time samples before and including the PC1 heating profile exceeding the threshold
%         P1pm = squeeze(mean(P1(:,:),2)); % pm = pressure mean
%         pause
%      end
%   end
%end % basin index for loop 
% Calc the spectrum averaged over i-lag time index, lag the big PC1 index
%for i=1:nlag+1 % time mean over i=1=nlag multiple indices
%   P1_ilag_tm(:,i) = mean(P1pm(:,i:nlag+1:end),2); % tm = time mean
%end
   %x = [0:dx:(nx-1)*dx]'; 
%xx = repmat(x,1,nz); % repmat(x,n,m) repeats matrix x in row (n-times) and column (m-times)
%X = repmat(x,1,nx);
%Z = linspace(0,max(z(:,nz)),nx); % resample z at the smallest height, which is 66 meters, hence about 3200 grid points vertically
%Z = repmat(Z,nx,1);
%ke_intrp=interp2(xx,z,ke,X,Z);
%contourf(1:nx,crmlev,ke')
%set(gca,'YDir','Reverse')
% FFT for Kinetic Energy
%KE = fftn(ke); % a n-dimensional fft, in my case, just doing fft in x-dir and then in z-dir
%
%P1 = P2(1:nx/2+1,1:nz/2+1); % k = 0 to nx/2 and 0 to nz/2 is needed, the other half has symmetric amplitude but flipped left-right.
%P1(2:end-1,:) = 2/nx*P1(2:end-1,:); % the scaling factors for waves k = 1 to nx/2-1 are all 2/nx
%P1(1,:)       = 1/nx*P1(1,:); % the k = 0 and nx/2 corresponding to cos(0) and cos(2pi) waves both have scaling factor of 1/nx
%P1(end,:)     = 1/nx*P1(end,:);
%P1(:,2:end-1) = 2/nz*P1(:,2:end-1);
%P1(:,1)       = 1/nz*P1(:,1);
%P1(:,end)     = 1/nz*P1(:,end);
%diag(P1); % select the frequencies k_zonal = k_vertical (k_vertical is per
%contourf(P1')
%
%


% CREATE INDEX SET TO CRM_U W
% CREATE A LS of filenames and go through line by index number
%CHECK THE PRESSURE LEVELS OF CRM



