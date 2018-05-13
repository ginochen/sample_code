function plot_mcs_size_dist(season,dosave)
%function plot_mcs_size_dist(season)
% season={'DJF','JJA','SON','MAM'}; for is=1:4; plot_mcs_size_dist(season{is},dosave); end
% size is plotted in scatter plot with darker color representing larger size, max size is around 16 points
%dosave = 1;
docbar=0;
docdf=0;
figure('units','normalized','outerposition',[0 0 1 1])
res = 'f09f09';
%casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch';
if ismember(season,{'MAM','SON'})
  casei='F_2000_SPCAM_m2005_3hrly1';
else
  casei='F_2000_SPCAM_m2005_3hrly2';
end
% new diro
dlat = 'lat9090'
load(['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season '/mcs_clusters_1.mat']);
diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/'];
% old diro
%load(['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/' res '/' season '/mcs_clusters.mat']);
%diro = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season '/'];
dircbar = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/'];
sdeg = 0.9414*1.25;
%dlon = diff(lon(1:2));
%dlat = diff(lat(1:2));
%area = dlon;
%rad = 6378e3;
%cir = 2*pi*rad;
mcsnllmax = max(mcsnll,[],2);
ns = max(mcsnllmax);
%ns = round(max(mcsnllave)) % maximum cluster size
%ns4plot = ns 
ns4plot = 17
cmapp = cmap(18); % 180 rows
icmap = round(linspace(1,size(cmapp,1),3))
icmap(2) = 4;
coast_centered(0);
for is = 1:ns % # of points
  sCl{is} = find(mcsnllmax==is); % cluster indices assoc to is-(# of points) 
%  sCl{is} = find(abs(mcsnllave-is)<0.5); % cluster indices assoc to is-(# of points) 
  if isempty(sCl{is}); continue; end
  for i = 1:numel(sCl{is})
    tmp{is}(i,1:2) = mcsilltcentroids{sCl{is}(i)}(1,1:2); % initial ilon-ilat of sCl{is}-cluster
  end
%  plot(lon(tmp{is}(1:2,1)),lat(tmp{is}(1:2,2)),'o','linewidth',2);pause
  if is<=2
    iss=1;
%  elseif is<=2
%    iss=2;
  elseif is<=5
    iss=2;
  else
    iss=3;
  end
  scatter(lon(tmp{is}(:,1)),lat(tmp{is}(:,2)),(is^1.3)*5,'o','linewidth',2,'MarkerEdgeColor',cmapp(icmap(iss),1:3)); 
  %plot(lon(tmp{is}(:,1)),lat(tmp{is}(:,2)),'Markersize',is*5,'o','linewidth',2,'MarkerEdgeColor',cmapp(icmap(is),1:3)); 
  hold on; 
end
%ylim(latlim);
ylim([-70 70])
xlim([min(lon),max(lon)]);
xticks([0:50:350])
set(gca,'xticklabel',{'0E','50E','100E','150E','160W','110W','60W','10W'})
%h = colorbar('Ticks',[5,9,13,17],'TickLabels',{'1','2','5','17'});
h = colorbar('Ticks',[6.3,11.6,17],'TickLabels',round(sdeg*[2,5,17]));
colormap(cmapp(icmap,1:3));
caxis([1 ns4plot]); % 3 hour is the minimum, and 22*3 is the maximum hours
set(get(h,'title'),'string','[deg^2]');
set(gca,'color','w')
set(gca,'fontsize',30)  
set(gca,'box','on')
colorbar('off')
daspect([.30 .5 1])
%title(['MCS Cluster Time-Averaged Size for ' season])
%title([season ' (n = ' num2str(nCl) ')'])
title([season])
if dosave
  disp(diro)
  saveas(gcf,[ diro '/cluster_tave-size.fig'])
  saveas(gcf,[ diro '/cluster_tave-size.png'])
  crop([diro '/cluster_tave-size.png'])
else
  sprintf(['saveas(gcf,' diro '/cluster_tave-size.fig'])
  sprintf(['saveas(gcf,' diro '/cluster_tave-size.png'])
end
disp('saved!')
if docbar
figure('units','normalized','outerposition',[0 0 1 1])
axis off
h = colorbar('southoutside','Ticks',[6.3,11.6,17],'TickLabels',round(sdeg*[2,5,17]));
set(gca,'Fontsize',44); 
colormap(cmapp(icmap,1:3));
caxis([1 ns4plot]); % 3 hour is the minimum, and 22*3 is the maximum hours
h.Label.String = '[deg^2]';
h.Position = [0.05 0.1170 0.7750 0.0181];
h.Label.Position = [19 0 0];
if dosave
  saveas(gca,[dircbar '/cluster_tave-size_colorbar.fig'])
  saveas(gca,[dircbar '/cluster_tave-size_colorbar.png'])
  crop([dircbar '/cluster_tave-size_colorbar.png'])
end
end
%return
%
if docdf
figure('units','normalized','outerposition',[0 0 1 1])
[f x]= ecdf(mcsnllave);
plot(x*sdeg,f*100,'o','MarkerEdgeColor','b','MarkerFaceColor','b','Markersize',10); yticks(0:20:100); % Empirical CDF of the tave-size in units of cell 
xlim([0,16])
%set(gca,'XtickLabel',[0:2:16])
set(gca,'fontsize',44)  
%title(['MCS Cluster Time-Averaged Size CDF for ' season])
%title([season ' (n = ' num2str(nCl) ')'])
xlabel('mean size [deg^2]')
ylabel('CDF [%]')
set(gcf,'color','w')
if dosave 
  saveas(gcf,[ diro '/CDF_tave-size.fig'])
  saveas(gcf,[ diro '/CDF_tave-size.png'])
  crop([diro '/CDF_tave-size.png'])
else
  sprintf(['saveas(gcf,' diro '/CDF_tave-size.fig'])
  sprintf(['saveas(gcf,' diro '/CDF_tave-size.png'])
end
end
%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%
% season = {'MAM','JJA','SON','DJF'}; for is=1:4; plot_mcs_size_dist(season{is}); close all; end
%
%
%
%

%plot(mcsilltcentroids{(sCl{1})
%for ic = 1:nCl
%  tmp(ic,1:3) = [lon(mcsilltcentroids{ic}(1,1)), lat(mcsilltcentroids{ic}(1,2)),mcsnllave(ic)];
%end
%latt = lat(lat>=-25 & lat <=25);
%[LON LAT] = meshgrid(lon,lat);
%Z = griddata(tmp(:,1),tmp(:,2),tmp(:,3),LON,LAT,'nearest');
%mesh(LON, LAT, Z); axis tight; 
%
%surf(tmp(:,1),tmp(:,2),tmp(:,3));
%mCl = find(mcsnllave >3 & mcsnllave <=6);
%bCl = find(mcsnllave >6);
%for i = 1:numel(mCl)
%  tmpm(i,1:2) = mcsilltcentroids{mCl(i)}(1,1:2);
%end
%for i = 1:numel(bCl)
%  tmpb(i,1:2) = mcsilltcentroids{bCl(i)}(1,1:2);
%end
%for i = 1
%scatter(tmps(:,1),tmps(:,2),'k'); hold on;
%scatter(tmpm(:,1),tmpm(:,2),'b'); hold on;
%scatter(tmpb(:,1),tmpb(:,2),'r'); hold on;
