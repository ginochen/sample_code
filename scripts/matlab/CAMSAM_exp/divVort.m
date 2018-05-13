function divVort(it,fname)
tic
run /nethome/gchen/scripts/matlab/startup.m
fname
%%%%%%%%%%%%%%%%%%%%%%%%
%    Job conditions    %
%%%%%%%%%%%%%%%%%%%%%%%%
dofft=0;
dodivvort=0;
%%%%%%%%%%%%%%%%%%%%%
%      Paths        %
%%%%%%%%%%%%%%%%%%%%%
% CAUTION!!! spcam_actual_m2005_f09f09_branch.cam.r.*.nc starts from
% 0001-01-14-01800 instead of 0001-01-14-03600 which the rh0 variables
% starts from
Case.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
%Case.dir    = '/projects/rsmas/kirtman/gchen/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
dir_r       = dir([Case.dir '/spcam_actual_m2005_f09f09_branch.cam.r.*-*-*-*.nc']);
Case.name   = {dir_r.name};
archive        = ['/projects/rsmas/kirtman/gchen/archive/'];
camCase.dir    = [archive 'F_2000_4SPCAM_m200501/atm/rest/'];
spcamCase.dir  = [archive 'spcam_cam_rh0_m2005/atm/rest/'];
camCase.name   = 'F_2000_4SPCAM_m200501.cam.rh0';
spcamCase.name = 'spcam_actual_m2005_f09f09_branch.cam.rh0';
landfrac  = loadvar('LANDFRAC',camCase,1,0);
%%%%%%%%%%%%%%%%%%%%%%%%
%    Load dimensions   %
%%%%%%%%%%%%%%%%%%%%%%%%
CaseDim.name   = 'spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-01800';
CaseDim.dir    = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/spcam_actual_m2005_f09f09_branch/atm/rest/';
lat  = loaddim('lat',CaseDim);
lon  = loaddim('lon',CaseDim);
lev  = loaddim('lev',CaseDim);  % level centers
ilev = loaddim('ilev',CaseDim); % level interfaces, two interfaces surround one level center
                                % index 1 (TOA pressure) to nlev (surface pressure) 
%%%%%%%%%%%%%%%%%%%%%%%%
%     Parameters       %
%%%%%%%%%%%%%%%%%%%%%%%%
dlev = ilev(2:end)-ilev(1:end-1); % level thickness (>0)
nlat=numel(lat);
nlon=numel(lon);
nlev=numel(lev);
nx    = 32; % evenly spaced
nz    = 28; % levels of CAM from surface to nz levels (according to cam/crm_physics.F90 line 992)
km2m  = 1000; % 1km = 1000m
dx    = 4*km2m; % 4 km grid box
%nlag  = 7; % extract the nlag timesteps for each grid point that exceeds the PC1 threshold at 0-lag timestep
nlag = 48;
lagind = 0:-10:-nlag;%nlag:-1:-nlag; % actual lag index 
%lagind_extra = nlag+1:-1:-nlag; % extra lag is to get the instantaneous field by stracting previous timestep
% dtlag = 30; % 1800 sec = 30 min for CAM timestep
% timelag = nlag*dtlag:-dtlag:0; % 0 lag is the big PC1 event timestep
g = 9.8; 
crmlev   = lev(nlev:-1:nlev-nz+1)*1e2; % [Pa] 1 (surface pressure) to nlev (TOA)
%crmdp = ilev(2:nz+1) - crmilev(1:nz);
crmdlev  = dlev(nlev:-1:nlev-nz+1)*1e2;
%crmdp(1) = crmlev(1)-ilev(end); % p-diff between 1st level to bottom-most pressure
%crmdp(2:nz)= crmlev(2:nz) - crmlev(1:nz-1); % the center to center pressure thickness of each level
%crmdlev = dlev(nlev:-1:nlev-nz+1);
%
%zint = linspace(1000,31000,16)'; % interpolated height with 2km grid-spacing
nzi  = 10; % # of interpolated z levels, nzi=61(~500m), nzi=180(~150m)
maxz = 18000; % ~31500m are most of the heighest heights of each vert columns
minz = 200; % ~ 55m are most of the lowest heights of each vert columns
zint = linspace(minz,maxz,nzi)'; % interpolated height with ~500m grid-spacing
dzi  = zint(2)-zint(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Heating profile SVD 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
threshold_bigPC = 10; % the 1st EOF mode largest magnitude vertically is 0.5, and the average PC has a vertical maximum of 5 K/day heating rate, so a threshold of 10 to get 10*0.5=5 is necessary
iEOF = 1; % EOF1 index
nzz=17;
start = [1 1 1]; % starting index for lon(1:288) lat(1:192) nx*nz(1:896) (maybe use lat(76:117) 20 degrees north-south?) 
count = [inf inf inf]; % number of index to read in var, inf means use all elements in that dimension
bsign = [1 -1 1 1]; % the sign for PC1 in different basin since basin 2 has a different signed distribution
itlag = 1; % counter for lag index
for ilag=lagind % lagging big PC1's index by ilag, so ilag=0 is the big PC1 time index, negative values are lead
   u  = load_reshape(Case, it-ilag, 'CRM_U',  start, count, nlon, nlat, nx, nz);
   w  = load_reshape(Case, it-ilag, 'CRM_W',  start, count, nlon, nlat, nx, nz);
   T  = load_reshape(Case, it-ilag, 'CRM_T',  start, count, nlon, nlat, nx, nz);
   qT = load_reshape(Case, it-ilag, 'CRM_QT', start, count, nlon, nlat, nx, nz);
   qc = load_reshape(Case, it-ilag, 'CRM_QC', start, count, nlon, nlat, nx, nz);
   qi = load_reshape(Case, it-ilag, 'CRM_QI', start, count, nlon, nlat, nx, nz);
   qr = load_reshape(Case, it-ilag, 'CRM_QR', start, count, nlon, nlat, nx, nz);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %  Large PC1 lon-lat index selection  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if (ilag == lagind(1)) % nlag is the first loop for ilag
      for ib = 1:4 % basin index 
         iitmp = iillt_zm{ib}(:,find(bsign(ib)*pc_spz{ib}(:,iEOF) > threshold_bigPC )); % pick out the big PC1 amplitude lon, lat, time indices and use it to study the CRM KE spectrum
                                                                                        % rows: 1(lon) 2(lat) 3(time)
         iillt_zm_PC{ib} = sortrows(iitmp',3)'; % sort the rows according to the 3rd (time) column of iitmp'<--- this sorting might not be necessary
         iCAMll{ib} = find(iillt_zm_PC{ib}(3,:)==it); % find the column indices of iillt_zm_PC{ib} for lat & lon at it-timestep, use these lat & lon for all lag-lead timesteps
      end
   end
   for ib = 1:4
      if (isempty(iCAMll{ib})==0) % if not empty 
         ill_CAM = 1; % counter for the lon lat indices for CAM
         for ill = iCAMll{ib}; % lon lat index 
            ilon = iillt_zm_PC{ib}(1,ill);
            ilat = iillt_zm_PC{ib}(2,ill);
            ilons = getlongs(ilon,nlon,ilat,landfrac); 
            if ( ilons ~= 0); % if there are land covered, discard this lon-lat index
               iilons=1; %ilons counter
               for ilon = ilons
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %      interpolate variable (p to z)     %
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  for ix = 1:nx
                     for iz = 1:nz
                        [rho(iz) qv(iz)] = density_temp(T(ilon,ilat,ix,iz),crmlev(iz),qT(ilon,ilat,ix,iz),qc(ilon,ilat,ix,iz),qi(ilon,ilat,ix,iz));
                        dz(iz)  = crmdlev(iz)/(rho(iz)*g); %dp = rho * g * dz
                        z(iz) = sum(dz(1:iz));
                     end;
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %             buoyancy freq                   %
                     % (bring a parcel from lower to higher level) %
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     n2 = - g^2 ./ rho(1:nz-1) .* (rho(2:nz)-rho(1:nz-1))./(z(2:nz)-z(1:nz-1)); % assume n2 misses the highest level, so it'll be at z(1) to z(nz-1)
%                     n2crm{ib}(ix,1:nzi,itlag,ill_CAM) = interp1qr(z(1:nz-1)', n2', zint);
                     var{ib}(iilons,ix,1:nzi,itlag,ill_CAM,:) = interp1qr(z',[squeeze(u(ilon,ilat,ix,1:nz)), ...
                                                                       squeeze(w(ilon,ilat,ix,1:nz)), ...
                                                                       squeeze(T(ilon,ilat,ix,1:nz)), ...
                                                                       squeeze(qT(ilon,ilat,ix,1:nz)), qv', ...
                                                                       squeeze(qr(ilon,ilat,ix,1:nz)), rho'], zint);
                  end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %  specific mass divergence  %
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  if (dodivvort)
                     [ div{ib}(1:nx,1:nzi,itlag,ill_CAM), vort{ib}(1:nx,1:nzi,itlag,ill_CAM) ] = div_vort(var{ib}(:,:,itlag,ill_CAM,1),var{ib}(:,:,itlag,ill_CAM,2),dx,dzi);
                  end         
                  %%%%%%%%%%%%%%%%%
                  %  ke spectrum  %
                  %%%%%%%%%%%%%%%%%
                  if (dofft)
                     FFTke{ib}(1:nx/2+1,1:nzi/2+1,itlag,ill_CAM) = fftke(var{ib}(:,:,itlag,ill_CAM,1),var{ib}(:,:,itlag,ill_CAM,2),nx,nzi);
                  end
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % rainfall gridpoint counts %
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %qr{ib}(itlag,ill_CAM) = numel(nonzeros(qr(ilon,ilat,:,1)))
                  %qrcam{ib}(iilons,1:nx,itlag,ill_CAM) = qr(ilon,ilat,:,1); %[kg/kg]
                  iilons=iilons+1;
               end
               ill_CAM = ill_CAM+1; % counter for the lon-lat CAM indices at a fixed time
            end
         end
      end
   end
   %contourf(reshape(permute(squeeze(qrcam{ib}(:,:,itlag,ill_CAM,:)),[1,3,2]),nx*(iilons-1),nz)','linestyle','none'); pause
   itlag = itlag+1;
end


%%%%%%%%%%%%%%%%%%%%%%%
%   Save variables    %
%%%%%%%%%%%%%%%%%%%%%%%
parm.var = {'u','w','T','qT','qv','qr','rho'}; % qc+qi = qT-qv 
parm.nx = nx; parm.zint = zint; parm.dzi=dzi; 
parm.nzi=nzi; parm.crmlev=crmlev; parm.nlag=nlag;
parm.lagind=lagind; parm.dx=dx; parm.iillt_zm_PC=iillt_zm_PC;
save(fname,'parm','var')
%save(fname,'parm','var','n2crm','div','vort','FFTke','qr')
%save(['T_qc_qt_qr_PC1_' num2str(it) '.mat'],'parm','Tcrm','qccrm','qTcrm','qrcrm')

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          FUNCTIONS USED BY divVort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [ilons]     = getlongs(ilon, nlon, ilat, landfrac)
   ilons = ilon-5:ilon+5;
   if (ilons(1)<=0);        ilons(find(ilons<=0))      = [nlon-sum(ilons<=0)+1:nlon]; end
   if (ilons(end)>=nlon+1); ilons(find(ilons>=nlon+1)) = [1:sum(ilons>=nlon+1)];      end
   if (sum(landfrac(ilons,ilat))~=0) % if there are land indices involved, discard this sample
      ilons = 0;
   end
