function latC = mcs_findlatzone(fin,season)
% find whether the clusters belong to tropics (0), midlat (1), or polar zone (2)
%casei='F_2000_SPCAM_m2005_3hrly2';
%diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season '/']; disp(diri)
%load([diri '/mcs_clusters_1.mat'],'lat','mcsilltcentroids','nt4Cl','nCl');
load(fin,'lat','mcsilltcentroids','nt4Cl','nCl','island');
%sum(island==0)
tratiothres = 0.9; % ratio of lifetime in tropics or subtropics
for ic=1:nCl
  lats = lat(mcsilltcentroids{ic}(:,2));
  if ( sum(lats>=-23.5 & lats<=23.5)/nt4Cl(ic)>=tratiothres ) % trops
    latC(ic) = 0;
  elseif ( sum(lats<=-23.5 & lats>=-66.5 | lats>=23.5 & lats<=66.5)/nt4Cl(ic)>=tratiothres ) % midlat
    latC(ic) = 1;
  else % polar or in-between trop & subtrop or subtrop & polar
    latC(ic) = 2;
  end
end
