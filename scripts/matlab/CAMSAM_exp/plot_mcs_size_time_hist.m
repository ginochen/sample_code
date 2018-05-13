function plot_mcs_size_time_hist(season)
% Fiolleau and Roca 2013 Figure 4
figure('units','normalized','outerposition',[0 0 1 1])
res = 'f09f09';
load(['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/' res '/' season '/mcs_clusters.mat']);
diro = ['/Users/g/archive/matlab/figure/mcs_cluster/' res '/' season]
%load('~/scripts/matlab/cmaps','cmap_mcs_size_time0')
sdeg = 0.9414*1.25;
xedges = [5:5:200];
yedges = [1:20];
histogram2((nt4Cl+2)*3,mcsnllave*sdeg,'XBinEdges',xedges,'YBinEdges',yedges,'DisplayStyle','tile','Normalization','probability')%,'binmethod','integers');
%scatter((nt4Cl+2)*3,mcsnllave*sdeg)
set(gca,'yscale','log','xscale','log')
corr(nt4Cl',mcsnllave')
pause
%colormap(cmap(19)); 
%colormap(cmap_mcs_size_time); 
%h= colorbar; 
%set(get(h,'title'),'string','[Cells]');
caxis([0 0.7])
%set(h,'Ticks',[0:400:6000]); 
xlim([3 200])
ylim([1 20]); 
set(gca,'Fontsize',40);
title(season)
xlabel('duration [h]')
ylabel('mean size [deg^2]')
saveas(gcf,[diro '/mcs_population_size_time.fig'])
saveas(gcf,[diro '/mcs_population_size_time.png'])
crop([diro '/mcs_population_size_time.png'])

%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%
% season = {'MAM','JJA','SON','DJF'}; for is=1:4; plot_mcs_size_time_hist(season{is}); close all; end
