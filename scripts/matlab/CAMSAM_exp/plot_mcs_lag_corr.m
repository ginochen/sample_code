function plot_mcs_lag_corr(season,ilndocn,dosave)
% season='JJA',ilndocn=0; dosave=0; plot_mcs_lag_corr(season,ilndocn,dosave)
docbar=1; % save colorbar?
%%{
%vn = {'LImin','dU03km_cam'}; % weak corr over ocean
%vn = {'dU03km_cam','LIave'}; % weak corr 
%vn = {'LImin','du03km_abs'}; % bad corr
%vname = {'shear','LI'};
%vn = {'dU03km_cam','precl'};
%vname = {'shear','stratiform rainrate'};
%vn = {'precc','precc'};
%vname = {'convective rainrate','convective rainrate'};
%vn = {'prect','dU03km_cam'};
%vname = {'Total rainrate','shear'};
%vn = {'du03km_abs','precc'}; % bad corr
%vn = {'dU03km_cam','precc'}; % slight pos corr over JJA ocn & lnd (shear leads to rain), bad in DJF 
%vname = {'shear','convective rainrate'};
%vn = {'dU03km_stormrel','precc'}; % slight pos corr over JJA ocn & lnd (shear leads to rain), bad in DJF 
%vname = {'storm shear','convective rainrate'};
%vn = {'dU03km_cam','Bmax'}; % lnd JJA better pos corr
%vname = {'shear','CPI'};
%vn = {'LIave','frac_precc'};
%vname = {'Convective rainrate fraction','CPI'};
%vn = {'prect','Bave'}; % no lag positively corr up to 0.5
%vn = {'prect','Bmax'}; % no lag positively corr up to 0.5
%vname = {'total rainrate','CPI'};
%vn = {'precc','Bave'}; % no lag positively corr up to 0.5
%vn = {'precc','Bmax'}; % no lag positively corr up to 0.5
%vname = {'convective rainrate','CPI'};
%vn = {'precl','Bmax'}; % no lag positively corr up to 0.5
%vname = {'stratiform rainrate','CPI'};
%vn = {'precl','precc'}; % lag positively since precl devel after precc
%vname = {'stratiform rainrate','convective rainrate'};
%vn = {'prect','precc'}; % lag positively since precl devel after precc
%vname = {'total rainrate','convective rainrate'};
%vn = {'dU03km_cam','prect'}; % no lag positively corr
%vname = {'du','Total rainrate'};
%vn = {'LIave','prect'};
%vname = {'LI','total rainrate'};
%vn = {'LIave','precc'};
%vn = {'LImin','precc'};
%vname = {'LI','convective rainrate'};
%vn = {'LIave','frac_precc'};
%vname = {'LI','Convective rainrate fraction'};
%%}
%%{
%vn = {'LIave','Bave'};
%vn = {'LImin','Bmax'};
%vn = {'LIave','Bmax'};
%vname = {'LI','CPI'};
%%}
%{
%vn = {'Bmax','du03km'};
%vname = {'B','shear'};
%}
%vn = {'LI','Bave'};
%vname = {'Lifted Index','B'};
lndocnstr = {'ocn','lnd'};
casei = 'F_2000_SPCAM_m2005_3hrly2';
cd(['/Users/g/archive/matlab/' casei '/atm/hist/lat2525/' season '/mcs_cluster_stats']);
diro = ['/Users/g/archive/matlab/' casei '/figure/lat2525/' season '/' lndocnstr{ilndocn+1} '/'];
%cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/DJF/mcs_cluster_stats
eval(sprintf(['load /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/' season '/mcs_cluster_stats/mcs_stats_' lndocnstr{ilndocn+1} ' nlndocnCl iclndocn '  ]))
eval(sprintf(['load /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/' season '/mcs_clusters']));
load(['mcs_stats_' lndocnstr{ilndocn+1} '.mat'],vn{1},vn{2});
itime = 8; %threshold time for cluster index
ii=1;
maxnt = 20
for nt2=3:maxnt
  if nlndocnCl{nt2} <= 1; continue; end % cannot calc stats with one sample
  for iic = 1:nlndocnCl{nt2}
    icc(ii) = iclndocn{nt2}(iic); % cluster index after removing end-time cluster indices
    ii=ii+1;
  end
end
tmp = (nt4Cl(icc)+2)*3;
idx = tmp>itime; % cluster index for long-lived clusters

eval(sprintf(['v1 = ' vn{1} '.mcscomposite.data(:,idx)'';']));
eval(sprintf(['v2 = ' vn{2} '.mcscomposite.data(:,idx)'';']));
if ismember(vn{1},{'du03km_cam_abs','du03km_abs','dU03km_cam','du03km'}) & ismember(vn{2},{'LImin','LIave'})
  if ilndocn==1
  ctick = [-0.15:0.1:0];
  climit = [-0.15 0];
  cm = flipud(cmap(19));
  else
  ctick = [-0.1:0.1:0.1];
  climit = [-0.1 0.1];
  cm = flipud(cmap(19));
   
  end
elseif ismember(vn{1},{'du03km_cam_abs','du03km_abs','dU03km_cam','du03km'}) & ismember(vn{2},{'Bmax','Bave'})
  ctick = [0:0.1:0.2];
  climit = [0 0.2];
  cm = cmap(19);
elseif ismember(vn{1},{'du03km_cam_abs','du03km_abs','dU03km_cam','du03km'}) & ismember(vn{2},{'precl'})
  ctick = [-0.1:0.1:0.1];
  climit =[-0.1 0.1];
  cm = cmap(19);
elseif ismember(vn{1},{'du03km_cam_abs','du03km_abs','dU03km_cam','du03km'}) & ismember(vn{2},{'precc'})
  ctick = [-0.1:0.1:0.1];
  climit =[-0.1 0.1];
  cm = cmap(19);
elseif ismember(vn{1},{'dU03km_stormrel'}) & ismember(vn{2},{'precc'})
  ctick = [-0.1:0.1:0.1];
  climit =[-0.1 0.1];
  cm = cmap(19);
elseif ismember(vn{1},{'precl'}) & ismember(vn{2},{'B','Bmax','Bave'}) 
  climit =[0.0 0.4];
  ctick = floor(1e2*linspace(0.0,0.4,3))*1e-2;
  cm = cmap(19);
elseif ismember(vn{1},{'prect'}) & ismember(vn{2},{'precc'}) 
  climit =[0.0 0.4];
  ctick = floor(1e2*linspace(0.0,0.4,3))*1e-2;
  cm = cmap(19);
elseif ismember(vn{1},{'precc'}) & ismember(vn{2},{'precc'}) % very similar to precc & Bmax
  ctick = [0:0.1:0.7]; 
  climit =[0.0 0.7];
  cm = (cmap(19));
elseif ismember(vn{1},{'precl'}) & ismember(vn{2},{'precc'}) 
  climit =[0.0 0.4];
  ctick = floor(1e2*linspace(0.0,0.4,3))*1e-2;
  cm = cmap(19);
elseif ismember(vn{1},'precc') & ismember(vn{2},{'Bave'}) 
  ctick = [0:0.1:0.5]; 
  climit =[0.0 0.5];
  cm = (cmap(19));
elseif ismember(vn{1},'precc') & ismember(vn{2},{'Bmax'}) 
  ctick = [0:0.1:0.7]; 
  climit =[0.0 0.7];
  cm = (cmap(19));
elseif ismember(vn{1},{'LIave','LImin'}) & ismember(vn{2},{'B','Bmax','Bave'}) 
  ctick = [-0.5:0.1:0];
  climit =[-0.5 0];
  cm = flipud(cmap(19));
elseif ismember(vn{1},{'LIave','LImin'}) & ismember(vn{2},{'precc'})
  climit =[-0.3 0.1];
  ctick = [-0.3:0.1:0.1];
%  ctick = floor(1e2*linspace(-0.3,0.1,3))*1e-2;
  cm = flipud(cmap(19));
end
figure('units','normalized','outerposition',[0 0 1 1])
[cc pv] = corr(v1,v2,'rows','complete');
lt =(linspace(0,100,20));
pv(pv<0.05)=1; % significant p values
pv(pv>0.05 & pv~=1)=0;
[iy ix] = find(pv);
[~,h]=contourf(lt,lt,cc,150,'linestyle','none'); hold on
scatter(lt(ix),lt(iy),'k*') % scatter plot points that are significant
caxis(climit);
plot(lt,[lt;lt]','k');
ylabel(['lifetime of ' vname{1} ' [%]'])
xlabel(['lifetime of ' vname{2} ' [%]'])
xticks(0:25:100);
xticklabels([0:25:100]);
yticks(0:25:100);
yticklabels([0:25:100]);
colormap(cm); hold on; 
set(gca,'Fontsize',44)
axis square
if dosave
  saveas(gcf,[diro '/corr_' vn{1} '_' vn{2} '.png'])
  crop([diro '/corr_' vn{1} '_' vn{2} '.png'])
end
disp(['saveas(gcf,' diro '/corr_' vn{1} '_' vn{2} '.png)'])

if docbar
  figure('units','normalized','outerposition',[0 0 1 1]); 
  axis off; 
  hc = colorbar('southoutside'); 
  set(hc,'ytick',ctick)
  colormap(cm); 
  caxis(climit);
 % colormap(cmap(ncmap(ivv(vname)),flipm(ivv(vname)),10)); 
 % caxis(caxislim{ivv(vname)}); 
  hc.Position = [0.05 0.2170 0.7750 0.0181];
  set(gca,'Fontsize',44); 
  hc.Label.String = 'corr'; 
  hc.Label.HorizontalAlignment = 'right';
  hc.Label.Position = [climit(2)+0.2*diff(climit) 0 0];
  saveas(gcf,['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/corr_' vn{1} '_' vn{2} '.png']); 
        crop(['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/corr_' vn{1} '_' vn{2} '.png'])
end
