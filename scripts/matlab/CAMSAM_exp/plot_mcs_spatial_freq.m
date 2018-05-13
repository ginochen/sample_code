function plot_mcs_spatial_freq(season,dosave)
% frequency = number of clusters per season per spatial point
%%%%%%%%%% MCS spatial frequency map %%%%%%%%%%%
% MAM JJA (summer) SON DJF (winter)

% go through each spatial points and count the time elements total in each season
% MAM '0001-03-01-0000' to '0001-05-31-75600'
if ~exist('dosave')
   dosave = 0
end
 figure('units','normalized','outerposition',[0 0 1 1])
% set( gca, 'Units', 'normalized', 'Position', [0.05 0 0.9 1] );
% axes('Units', 'normalized', 'Position', [0 0 1 1])
 res = 'f09f09';
 casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch'
% diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' season '/']
% diro = ['/Users/g/archive/matlab/' casei '/figure/' season '/'];
 diri = ['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/f09f09/' season '/'];
 diro = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season];
% dircbar = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/'];
 load([diri 'mcs_clusters.mat'],'latlim','nCl');
 if iscell(season) % climatology
   tmp = 0;
   for is = 1:4
     fo = [diri '/' season{is} '/mcs_spatialfreq.mat'];
     load(fo);
     tmp = tmp + mcs_spatialfreq;
   end
   mcs_spatialfreq = tmp;
 else
   fo = [diri  '/mcs_spatialfreq.mat'];
%   seasonTOD = containers.Map({'MAM','JJA','SON','DJF'},{{'0001-03-01-00000','0001-05-31-75600'},...
%         {'0001-06-01-00000','0001-08-31-75600'},{'0001-09-01-00000','0001-11-30-75600'},{'0001-12-01-00000','0002-02-28-75600'}});
   seasonTOD = containers.Map({'MAM','JJA','SON','DJF'},{{'0002-03-01-00000','0002-05-31-75600'},...
         {'0002-06-01-00000','0002-08-31-75600'},{'0002-09-01-00000','0002-11-30-75600'},{'0001-12-01-00000','0002-02-28-75600'}});
   is = seasonTOD(season);
 end
 condplot = 1; icmap = 12;
 if (~exist(fo))
   load([diri '/mcs_clusters.mat'],'mcsillt','t','lon','lat','latlim');
   [mcs_spatialfreq, LON, LAT] = mcs_spatial_freq(mcsillt,t,lon,lat,is{1},is{2});
   save(fo,'mcs_spatialfreq','LON','LAT')
   disp('rerun plot_mcs_spatial_freq.m to plot')
 else
   load(fo);
   ncon = 100;
   contourf(LON, LAT, mcs_spatialfreq, ncon, 'linestyle','none'); 
%   h = colorbar;
%   set(get(h,'title'),'string','[Cells]')
   cmapp = cmap(icmap);
   icmapp = round(linspace(1,size(cmapp,1),ncon));
   colormap(cmapp(icmapp,:));
   if iscell(season)>1
     caxis([0 30]); % 3 hour is the minimum, and 22*3 is the maximum hours
   else 
     caxis([0 10]);
   end     
   %set(get(h,'title'),'string','Frequency');
   hold on
   coast_centered(0); % center the coastline at 0 degrees
   ylim(latlim);
%   title([season])
   title([season ' (n = ' num2str(nCl) ')'])
   set(gca,'fontsize',30)  
%   set(gcf,'color','w')
   daspect([1 .5 1])
   if dosave
%     saveas(gcf,[diro '/mcs_population.fig'])
     saveas(gcf,[diro '/mcs_population.png'])
     crop([diro '/mcs_population.png'])
   end
disp('saved!')
pause
   figure('units','normalized','outerposition',[0 0 1 1])
   h = colorbar('southoutside'); 
   axis off
   h.Label.String = '[cells]';
   h.Position = [0.05 0.1170 0.7750 0.0181];
   h.Label.Position = [11 0 0];
   set(gca,'Fontsize',44);
%   set(get(h,'title'),'string','[Cells]')
   cmapp = cmap(icmap);
   icmapp = round(linspace(1,size(cmapp,1),ncon));
   colormap(cmapp(icmapp,:));
   if iscell(season)>1
     caxis([0 30]); % 3 hour is the minimum, and 22*3 is the maximum hours
   else 
     caxis([0 10]);
   end   
   dosave=0;
   display(['dosave for colorbar = ' num2str(dosave)])
   if dosave 
%     saveas(gca,[dircbar '/mcs_population_colorbar.fig'])
     saveas(gca,[dircbar '/mcs_population_colorbar.png'])
     crop([dircbar '/mcs_population_colorbar.png'])
   end
 end
% mcs_spatialfreqDJF = mcs_spatial_freq(mcsillt,t,lon,lat,'0001-12-01-00000','0001-02-28-75600');

disp([' dosave=0; season={''MAM'',''JJA'',''SON'',''DJF''}; for is=1:4; plot_mcs_spatial_freq(season{is},dosave); end']) 
