
% vertical velocity 
timesteps = 41:10:401;
ucam1D_w=[];
ucrm1D_w=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','varCRM2D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D','varCRM2D'},parm)
  end
  ucam1D_w = cat(2,ucam1D_w,varCAM1D(1:parm.nzi,:,idv.varCAM1D.w));
  ucrm1D_w = cat(2,ucrm1D_w,squeeze(mean(varCRM2D(1:parm.nx,1:parm.nzi,:,idv.varCRM2D.w),1)));
end
zint = parm.zint/1000; % m to km
ucam1D_w(ucam1D_w==0)=NaN;
ucrm1D_w(ucrm1D_w==0)=NaN;
ucam1D_w=ucam1D_w*1000; % from kg/kg to g/kg
ucrm1D_w=ucrm1D_w*1000;
maxw = 0.08; 
%maxw = 1.2; 
minw = 0.001;
nbins = 80
wcam_bincenters=[];
wcrm_bincenters=[];
npoint_wcam=[];
npoint_wcrm=[];
wcam_binedges = linspace(minw,maxw,nbins);
wcrm_binedges = linspace(minw,maxw,nbins);
figure
for i=1:nbins-1; wcam_bincenters(i) = mean(wcam_binedges(i:i+1));end;
for i=1:nbins-1; wcrm_bincenters(i) = mean(wcrm_binedges(i:i+1));end;
for iz=1:parm.nzi;
  xx = histogram(ucam1D_w(iz,:),wcam_binedges);
  npoint_wcam(1:nbins-1,iz) = xx.Values;
  xx = histogram(ucrm1D_w(iz,:),wcrm_binedges);
  npoint_wcrm(1:nbins-1,iz) = xx.Values;
end
Npoints_wcam = sum(npoint_wcam(:));
Npoints_wcrm = sum(npoint_wcrm(:));
red=[252, 187, 184
     252, 106, 98
     147, 7, 0]/255;
blue=[184, 196, 252
      98, 131, 252
      0, 24, 147]/255;
clev_w=[0.001 0.02 0.03];
linewidths=[2 3 5];
for i=1:numel(clev_w)
contour(wcrm_bincenters,zint, npoint_wcrm'/Npoints_wcrm,[clev_w(i) clev_w(i)],'color',blue(i,:),'linewidth',linewidths(i)); hold on
contour(wcam_bincenters,zint, npoint_wcam'/Npoints_wcam,[clev_w(i) clev_w(i)],'color',red(i,:),'linewidth',linewidths(i)); hold on
end
set(gcf,'color','w')
set(gca,'box','off')
set(gca,'Fontsize',30)
xlim([minw maxw])
xticks([0.001 0.01 0.02 0.03 0.04 0.05 0.06])
