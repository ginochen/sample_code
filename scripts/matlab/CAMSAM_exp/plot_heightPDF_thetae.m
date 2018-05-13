% thetae
% looks identical, needs to plot ice-\theta_e 
timesteps = 41:10:401;
ucam1D_thetae=[];
ucrm1D_thetae=[];
ucrm1D_qr=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D'},parm)
  end
  ucam1D_thetae = cat(2,ucam1D_thetae,varCAM1D(1:parm.nzi,:,idv.varCAM1D.thetae_cam));
  ucrm1D_thetae = cat(2,ucrm1D_thetae,varCAM1D(1:parm.nzi,:,idv.varCAM1D.thetae_crm));
%  ucrm1D_qr = cat(2,ucrm1D_qr,squeeze(mean(varCRM2D(1:parm.nx,1:parm.nzi,:,idv.varCRM2D.qr),1))); %include rain????????? maybe not since it peaks at no-cloud surface, and CAM thetae is in-cloud
end
zint = parm.zint/1000; % m to km
maxthetae = 410; 
minthetae = 310;
nbins = 80
thetaecam_bincenters=[];
thetaecrm_bincenters=[];
npoint_thetaecam=[];
npoint_thetaecrm=[];
thetaecam_binedges = linspace(minthetae,maxthetae,nbins);
thetaecrm_binedges = linspace(minthetae,maxthetae,nbins);
figure
for i=1:nbins-1; thetaecam_bincenters(i) = mean(thetaecam_binedges(i:i+1));end;
for i=1:nbins-1; thetaecrm_bincenters(i) = mean(thetaecrm_binedges(i:i+1));end;
for iz=1:parm.nzi;
  xx = histogram(ucam1D_thetae(iz,:),thetaecam_binedges);
  npoint_thetaecam(1:nbins-1,iz) = xx.Values;
  xx = histogram(ucrm1D_thetae(iz,:),thetaecrm_binedges);
  npoint_thetaecrm(1:nbins-1,iz) = xx.Values;
end
Npoints_thetaecam = sum(npoint_thetaecam(:));
Npoints_thetaecrm = sum(npoint_thetaecrm(:));
red=[252, 187, 184
     252, 106, 98
     147, 7, 0]/255;
blue=[184, 196, 252
      98, 131, 252
      0, 24, 147]/255;
clev_thetae=[0.001 0.02 0.03];
linewidths=[2 3 5];
for i=1:numel(clev_thetae)
contour(thetaecrm_bincenters,zint, npoint_thetaecrm'/Npoints_thetaecrm,10);hold on
contour(thetaecam_bincenters,zint, npoint_thetaecam'/Npoints_thetaecam,10);
%contour(thetaecrm_bincenters,zint, npoint_thetaecrm'/Npoints_thetaecrm,[clev_thetae(i) clev_thetae(i)],'color',blue(i,:),'linewidth',linewidths(i)); hold on
%contour(thetaecam_bincenters,zint, npoint_thetaecam'/Npoints_thetaecam,[clev_thetae(i) clev_thetae(i)],'color',red(i,:),'linewidth',linewidths(i)); hold on
end
set(gcf,'color','w')
set(gca,'box','off')
set(gca,'Fontsize',30)
xlim([minthetae maxthetae])
xticks([0.001 0.01 0.02 0.03 0.04 0.05 0.06])
