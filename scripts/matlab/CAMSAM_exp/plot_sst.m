function plot_sst(season)
filei = '/Users/g/draft/convection_stochastic/figure/sst_HadOIBl_bc_0.9x1.25_clim_c040926.nc'
sst = ncread(filei,'SST_cpl');
figure('units','normalized','outerposition',[0 0 1 1])
casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch';
load(['/Users/g/archive/matlab/' casei '/atm/hist/' season '/mcs_clusters.mat']);
diro = ['/Users/g/archive/matlab/' casei '/figure/' season '/'];
iseason = containers.Map({'JJA','DJF'},{[6:8],[1,2,12]});
lat(ilatlim)
contour(lon,lat(ilatlim(1):ilatlim(2)),median(sst(:,ilatlim(1):ilatlim(2),iseason(season)),3)',[26:1:30],'linewidth',2); hold on
coast_centered(0);
ylim(latlim);
xlim([min(lon),max(lon)]);
