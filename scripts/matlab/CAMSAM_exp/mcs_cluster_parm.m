function parm = mcs_cluster_parm(fin)
%function parm = mcs_cluster_parm(fin)
  p0   = ncread(fin,'P0'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  phis = ncread(fin,'PHIS'); % [Pa] in SI units not [hPa], can only use spcam_actual_m2005_f09f09_branch.cam.rh0.0001-01-14-03600, the other timesteps P0, PHIS are all zero
  hyam = flipud(ncread(fin,'hyam'));
  hybm = flipud(ncread(fin,'hybm')); 
  hyai = flipud(ncread(fin,'hyai'));
  hybi = flipud(ncread(fin,'hybi'));
  lat  = ncread(fin,'lat');
%  ilatlim = find(lat>=-20 & lat<=20)'; % tropical ilat
  lon  = ncread(fin,'lon');
  lev  = ncread(fin,'lev');  % level centers
  ilev = ncread(fin,'ilev'); % level interfaces, two interfaces surround one level center
  nz = 28;
  nx = 32;
  g  = 9.8; 
  km2m  = 1000; % 1km = 1000m
  dx    = 4*km2m; % 4 km grid box
  nzi  = 20; % # of interpolated z levels, nzi=61(~500m), nzi=180(~150m)
  maxz = 18000;
%  minz = 250; 
  minz = 100; 
  zint = linspace(minz,maxz,nzi)'; % interpolated height with ~500m grid-spacing
  dzi  = zint(2)-zint(1);
  i3km = find(abs(zint-3000) == min(abs(zint-3000))); 
  i6km = find(abs(zint-6000) == min(abs(zint-6000)));
  i10km = find(abs(zint-10000) == min(abs(zint-10000)));
  dz3km  = zint(i3km)  - zint(1); % the distance between first level to 3km, 3km-0km SHEAR defined in Jirak and Cotton
  dz6km  = zint(i6km)  - zint(1); % the distance between first level to 3km 
  dz10km = zint(i10km) - zint(1); % the distance between first level to 10km 
  parm.lon = lon; parm.lat = lat;
  R = 6.3781*10^6;
  ang = 7.2921159*10^-5;
  parm.lev = lev;
  parm.ilev = ilev;
  parm.dlon = 2*pi*R*cosd(parm.lat)*(lon(2)-lon(1))/360; % in meters
  parm.dlat = 2*pi*R*(lat(2)-lat(1))/360;
  parm.co = 2*ang*sind(parm.lat);
  parm.hyam = hyam; parm.hybm = hybm; parm.hyai = hyai; parm.hybi = hybi;
  parm.g = g; parm.nz = nz;
  parm.phis = phis; parm.p0 = p0; 
  parm.nx = nx; parm.zint = zint; parm.dzi=dzi;
  parm.i3km=i3km; parm.i6km=i6km; parm.i10km=i10km;
  parm.dz3km=dz3km; parm.dz6km=dz6km; parm.dz10km=dz10km;
  parm.nzi=nzi; 
  parm.dx=dx;
%  parm.season = season;
