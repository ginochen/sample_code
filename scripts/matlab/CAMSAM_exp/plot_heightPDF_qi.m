% qi
timesteps = 41:10:401;
ucam1D_qi=[];
ucrm1D_qi=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','varCRM2D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D','varCRM2D'},parm)
  end
  ucam1D_qi = cat(2,ucam1D_qi,varCAM1D(1:parm.nzi,:,idv.varCAM1D.icimr));
  ucrm1D_qi = cat(2,ucrm1D_qi,squeeze(mean(varCRM2D(1:parm.nx,1:parm.nzi,:,idv.varCRM2D.qi),1)));
end
zint = parm.zint/1000; % m to km
ucam1D_qi(ucam1D_qi==0)=NaN;
ucrm1D_qi(ucrm1D_qi==0)=NaN;
ucam1D_qi=ucam1D_qi*1000; % from g/g to g/kg
ucrm1D_qi=ucrm1D_qi*1000;
nbins=50;
maxqi = max(max(ucam1D_qi));
minqi = 1e-3;
%maxqi = 1e-02; 
%minqi = 2.6*1e-3;
qicam_bincenters=[];
qicrm_bincenters=[];
npoints_qicam=[];
npoints_qicrm=[];
qicam_binedges = linspace(minqi,maxqi,nbins);
qicrm_binedges = linspace(minqi,maxqi,nbins);
figure
for i=1:nbins-1; qicam_bincenters(i) = mean(qicam_binedges(i:i+1));end;
for i=1:nbins-1; qicrm_bincenters(i) = mean(qicrm_binedges(i:i+1));end;
for iz=1:parm.nzi;
  xx = histogram(ucam1D_qi(iz,:),qicam_binedges);
  npoints_qicam(1:nbins-1,iz) = xx.Values;
  xx = histogram(ucrm1D_qi(iz,:),qicrm_binedges);
  npoints_qicrm(1:nbins-1,iz) = xx.Values;
end
Npoints_qicam = sum(npoints_qicam(:));
Npoints_qicrm = sum(npoints_qicrm(:));
red=[252, 187, 184
     252, 106, 98
     147, 7, 0]/255;
blue=[184, 196, 252
      98, 131, 252
      0, 24, 147]/255;
clev_qi=[1e-3 2e-3 3e-3];
clev_qi=linspace(0.005,0.02,3);%clev_qi=linspace(0.005,0.025,3);%[0.005 0.015 0.02];
linewidths=[2 3 5];
for i=1:numel(clev_qi)
contour(qicam_bincenters,zint, npoints_qicam'/Npoints_qicam,[clev_qi(i) clev_qi(i)],'color',red(i,:),'linewidth',linewidths(i)); hold on
contour(qicrm_bincenters,zint, npoints_qicrm'/Npoints_qicrm,[clev_qi(i) clev_qi(i)],'color',blue(i,:),'linewidth',linewidths(i));
end
set(gcf,'color','w')
set(gca,'box','off')
set(gca,'Fontsize',30)
