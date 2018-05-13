% qc
timesteps = 41:10:401;
ucam1D_qc=[];
ucrm1D_qc=[];
ucrm1D_qr=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','varCRM2D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D','varCRM2D'},parm)
  end
  ucam1D_qc = cat(2,ucam1D_qc,varCAM1D(1:parm.nzi,:,idv.varCAM1D.icwmr));
  ucrm1D_qc = cat(2,ucrm1D_qc,squeeze(mean(varCRM2D(1:parm.nx,1:parm.nzi,:,idv.varCRM2D.qc),1)));
%  ucrm1D_qr = cat(2,ucrm1D_qr,squeeze(mean(varCRM2D(1:parm.nx,1:parm.nzi,:,idv.varCRM2D.qr),1))); %include rain????????? maybe not since it peaks at no-cloud surface, and CAM qc is in-cloud
end
%ucrm1D_qc = ucrm1D_qr + ucrm1D_qc;
zint = parm.zint/1000; % m to km
ucam1D_qc(ucam1D_qc==0)=NaN;
ucrm1D_qc(ucrm1D_qc==0)=NaN;
ucam1D_qc=ucam1D_qc*1000; % from kg/kg to g/kg
ucrm1D_qc=ucrm1D_qc*1000;
maxqc = 0.08; 
%maxqc = 1.2; 
minqc = 0.001;
nbins = 80
qccam_bincenters=[];
qccrm_bincenters=[];
npoint_qccam=[];
npoint_qccrm=[];
qccam_binedges = linspace(minqc,maxqc,nbins);
qccrm_binedges = linspace(minqc,maxqc,nbins);
figure
for i=1:nbins-1; qccam_bincenters(i) = mean(qccam_binedges(i:i+1));end;
for i=1:nbins-1; qccrm_bincenters(i) = mean(qccrm_binedges(i:i+1));end;
for iz=1:parm.nzi;
  xx = histogram(ucam1D_qc(iz,:),qccam_binedges);
  npoint_qccam(1:nbins-1,iz) = xx.Values;
  xx = histogram(ucrm1D_qc(iz,:),qccrm_binedges);
  npoint_qccrm(1:nbins-1,iz) = xx.Values;
end
Npoints_qccam = sum(npoint_qccam(:));
Npoints_qccrm = sum(npoint_qccrm(:));
red=[252, 187, 184
     252, 106, 98
     147, 7, 0]/255;
blue=[184, 196, 252
      98, 131, 252
      0, 24, 147]/255;
clev_qc=[0.001 0.02 0.03];
linewidths=[2 3 5];
for i=1:numel(clev_qc)
contour(qccrm_bincenters,zint, npoint_qccrm'/Npoints_qccrm,[clev_qc(i) clev_qc(i)],'color',blue(i,:),'linewidth',linewidths(i)); hold on
contour(qccam_bincenters,zint, npoint_qccam'/Npoints_qccam,[clev_qc(i) clev_qc(i)],'color',red(i,:),'linewidth',linewidths(i)); hold on
end
set(gcf,'color','w')
set(gca,'box','off')
set(gca,'Fontsize',30)
xlim([minqc maxqc])
xticks([0.001 0.01 0.02 0.03 0.04 0.05 0.06])
