function plot_mcs_time_dist(season,dosave)
%function plot_mcs_time_dist(season)
%dosave = 1;
docbar = 1;
docdf = 0;
figure('units','normalized','outerposition',[0 0 1 1])
res = 'f09f09';
thres_it=4; % threshold time length, 4=(4-2)*3=6hours
%casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch';
if ismember(season,{'MAM','SON'})
  casei='F_2000_SPCAM_m2005_3hrly1';
else
  casei='F_2000_SPCAM_m2005_3hrly2';
end
% new diro
dlat='lat9090';
load(['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season '/mcs_clusters_1.mat']);
diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/'];
% old diro
%load(['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/' res '/' season '/mcs_clusters.mat']);
%diro = ['/Users/g/archive/matlab/figure/mcs_cluster/' res '/' season];
dircbar = ['/Users/g/archive/matlab/figure/mcs_cluster/' res '/'];
nt = max(nt4Cl);
%nt = size(t4Cl,2); % t4Cl saves the unique time indices for each cluster, hence the columns of t4Cl shows the maximum time index length
doLong = 1; % plot the long lasting and short lasting MCS clusters only
%nt4plot = nt % corresponding to 51 hrs 
nt4plot = 60 % corresponding to 51 hrs 
cmapp = cmap(17); % 180 rows
if doLong
  icmap = round(linspace(1,size(cmapp,1),2))
else
  icmap = round(linspace(1,size(cmapp,1),nt4plot));
end
coast_centered(0);
for it = 1:nt % cluster lifetime
  if (numel(nt4ClRowInd{it})) % # of clusters that has lifetime it
    for i = 1:numel(nt4ClRowInd{it})
      tmp{it}(i,1:2) = mcsilltcentroids{nt4ClRowInd{it}(i)}(1,1:2); % starting centroid lon lat for cluster lifetime it
    end
    if doLong
      if it<=thres_it % length of cluster time
        itt=1;
%      elseif it<=12
%        itt=2;
      else
        itt=2;
      end
    else
      itt = it;
      if it > nt4plot % fix the maximum color to 66 hrs, so any cluster greater than 66 hrs will be plotted at the same color 
        itt = nt4plot; 
      end
    end
    scatter(lon(tmp{it}(:,1)),lat(tmp{it}(:,2)),(it^1.1)*5,'o','linewidth',2,'MarkerEdgeColor',cmapp(icmap(itt),1:3));%[1 1 1]*(nt-it)/nt);
    hold on;
  end
end
ylim([-70 70])
%ylim(latlim);
xlim([min(lon),max(lon)]);
xticks([0:50:350])
set(gca,'xticklabel',{'0E','50E','100E','150E','160W','110W','60W','10W'})
%h = colorbar('Ticks',[61,121,180],'TickLabels',{'10','30','174'}); % MCC lasts about 10 hours
%h = colorbar('Ticks',[61,121,180],'TickLabels',{num2str((thres_it-2)*3),'30','174'}); % MCC lasts about 10 hours
%colormap(cmapp(icmap,1:3));
%caxis([3 nt4plot*3]); % 3 hour is the minimum, and 22*3 is the maximum hours
%set(get(h,'title'),'string','[h]');
if ~docbar
  colorbar('off')
%caxis([3/24 nt4plot*3/24]); % 3 hour is the minimum, and 22*3 is the maximum hours
end
set(gcf,'color','w')
set(gca,'fontsize',30)  
set(gca,'box','on')
daspect([.30 .5 1])
%title([season ' (n = ' num2str(nCl) ')'])
title([season])
if dosave
  saveas(gcf,[diro '/cluster_lifetime.fig'])
  saveas(gcf,[diro '/cluster_lifetime.png'])
  crop([diro '/cluster_lifetime.png'])
else
  sprintf(['saveas(gcf,' diro '/cluster_lifetime.fig'])
  sprintf(['saveas(gcf,' diro '/cluster_lifetime.png'])
end
disp('saved!')
if docbar
figure('units','normalized','outerposition',[0 0 1 1])
axis off
h = colorbar('southoutside','Ticks',[91,180],'TickLabels',{num2str((thres_it-2)*3),num2str((nt-2)*3)}); % MCC lasts about 10 hours
colormap(cmapp(icmap,1:3));
set(gca,'Fontsize',44); 
h.Label.String = '[h]';
h.Label.Position = [195 0 0];
caxis([3 nt4plot*3]); % 3 hour is the minimum, and 22*3 is the maximum hours
if dosave  
  saveas(gca,[dircbar '/cluster_lifetime_colorbar.fig'])
  saveas(gca,[dircbar '/cluster_lifetime_colorbar.png'])
  crop([dircbar '/cluster_lifetime_colorbar.png'])
end
end
%
%return
if docdf
figure('units','normalized','outerposition',[0 0 1 1])
[f x]= ecdf(nt4Cl*3);
plot(x,f*100,'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',10);
yticks(0:20:100); % Empirical CDF of the liftime
set(gca,'fontsize',44)  
%title(['MCS Cluster Lifetime CDF for ' season])
%title([season ' (n = ' num2str(nCl) ')'])
%set(gca,'xticklabel',[0:20:200])
%xlim([0 200])
xlabel('duration [h]')
ylabel('CDF [%]')
%dosave = 1;
if dosave
  saveas(gcf,[diro '/CDF_lifetime.fig'])
  saveas(gcf,[diro '/CDF_lifetime.png'])
  crop([diro '/CDF_lifetime.png'])
else
  sprintf(['saveas(gcf,' diro '/CDF_lifetime.fig'])
  sprintf(['saveas(gcf,' diro '/CDF_lifetime.png'])
end
end
%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%
% season = {'MAM','JJA','SON','DJF'}; for is=1:4; plot_mcs_time_dist(season{is}); close all; end
