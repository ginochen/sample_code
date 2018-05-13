
%
spd = 86400; % seconds per day
% spcam m2005 case: 0001-01-14-03600 to 0001-01-29-63000
% cam m2005 case: 0001-01-14-03600 to 0001-01-29-64800
% Cam rh0 files has one additional timestep at the last step due to
% restarting at 63000 from spcam m2005 
%ivars_3d = [1,15,18,30,31,32];
%ivars_2d = [4,10,11,26];
%ivarsall = [ivars_3d,ivars_2d];
%ie = find(ivarsall==32); % index of e or e_abs
%ine= find(ivarsall~=32); % index of the others except e or e_abs
%
Case.dir = camCase.dir;
Case.name = [camCase.name '.0001-01-29-64800'];
slat = loaddim('slat',Case);
lat  = loaddim('lat',Case);
lon  = loaddim('lon',Case);
lev  = loaddim('lev',Case); % level centers
ilev = loaddim('ilev',Case); % level interfaces, two interfaces surround one level center
nlat=numel(lat);
nlon=numel(lon);
nlev=numel(lev);
dlev = ilev(2:end)-ilev(1:end-1); % level thickness
cutoff_lat = 20;
ilat = find(lat<=cutoff_lat & lat>=-cutoff_lat); % selected lat indices
ilon = [1:nlon]; % selected lon indices
%ilon = find(lon<=300 & lon>=120); % selected lat indices
itimes=[1:700]; % selected time indices
iilev = [1:nlev]; % selected level indices
%levels=[200,400,850,1000];
%levels = [0,120,350,820,930,957,976,992];
%for i=1:numel(levels); iilev(i) = find(abs(lev-levels(i)) == min(abs(lev-levels(i)))); end
%iilev(1)=0; % if levels(1)=0, then set iilev(1)=0
%nbin = 100;
%dtcondbin=linspace(-20,20,nbin);
