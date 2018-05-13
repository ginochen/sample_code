function plot_CF(season,ilndocn)
% contoured frequency for any mcs lifetime interpolated variables
% data(1:nt,1:nsamp)
%season = 'JJA';
%season = 'DJF';
%ilndocn = 1;
iv = 1;
vn = {'dU03km_cam','du03km_cam_abs','du03km_abs','LIave','Bave','Bmax','precc'};
%vn = {'LI','du03km_abs'};
vname = {'shear','shear_cam_u','shear_CRM', 'LI', 'B','Bmax','conv_prec'};
vrange = {[0 35],[0 35],[0 35],[-6 -2],[0 60],[0 300],[0 15]}; 
bwidth = [0.2,0.2,0.2,0.05,0.3,0.5,0.3];
lndocnstr = {'ocn','lnd'};
casei = 'F_2000_SPCAM_m2005_3hrly2';
cd(['/Users/g/archive/matlab/' casei '/atm/hist/lat2525/' season '/mcs_cluster_stats']);
diro = ['/Users/g/archive/matlab/' casei '/figure/lat2525/' season '/' lndocnstr{ilndocn+1} '/'];
%cd /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/DJF/mcs_cluster_stats

load(['mcs_stats_' lndocnstr{ilndocn+1} '.mat'],vn{iv});

eval(sprintf(['data = ' vn{iv} '.mcscomposite.data;']));
nanmean(data(:))
  nbin = diff(vrange{iv})/bwidth(iv);
  nt = size(data,1);
  hp = [1:nt];
  caxisrange =  [0 35];
  clev = linspace(caxisrange(1),caxisrange(2),100);
  [pdf, cdf, vbincenters] = binVarWithHeight(data,vrange{iv},nbin,nt);
  figure('units','normalized','outerposition',[0 0 1 1])
  hold on;
  contourf(vbincenters,hp,pdf'*100,clev,'linestyle','none'); 
  plot(nanmean(data,2),hp,'k');
  title([season ' ' lndocnstr{ilndocn+1}])
