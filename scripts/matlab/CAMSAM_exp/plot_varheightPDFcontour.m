function plot_varheightPDFcontour(varname,varmax,varmin)
%function plot_varheightPDFcontour(varname,varmax,varmin)
% for relative humidity inputs: 'relhum', 100, 0
figure
cd /Users/g/archive/matlab/figure/var/mcsVars/mcc % this is where mcc data are stored
timesteps = 41:10:401;
red=[252, 187, 184
     252, 106, 98
     147, 7, 0]/255;
blue=[184, 196, 252
      98, 131, 252
      0, 24, 147]/255;
nbins = 80;
clev=[0.001 0.003 0.006]; % contour levels
linewidths=[2 3 5];
% load varCAM1D into one single vector
[vcam, vcrm, nz, zint] = catVar(varname,timesteps);
% probability of variable with height
[prob_cam, vcam_bincenters] = bin_var_with_height(vcam,varmax,varmin,nbins,nz);
[prob_crm, vcrm_bincenters] = bin_var_with_height(vcrm,varmax,varmin,nbins,nz);
% plot var-height PDF contour
for i=1:numel(clev)
  contour(vcrm_bincenters,zint, prob_crm,[clev(i) clev(i)],'color',blue(i,:),'linewidth',linewidths(i)); hold on
  contour(vcam_bincenters,zint, prob_cam,[clev(i) clev(i)],'color',red(i,:),'linewidth',linewidths(i)); hold on
end
set(gcf,'color','w')
set(gca,'box','off')
set(gca,'Fontsize',30)
xlim([varmin varmax])
xticks([10:10:100])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prob_vcam, vcam_bincenters] = bin_var_with_height(vcam,varmax,varmin,nbins,nz)
% vcam(1:nz,1:nsamp)
% bin the variable at each height 
vcam_binedges = linspace(varmin,varmax,nbins);
for i=1:nbins-1; 
  vcam_bincenters(i) = mean(vcam_binedges(i:i+1));
end;
for iz=1:nz;
  xx = histogram(vcam(iz,:),vcam_binedges);
  npoint_vcam(1:nbins-1,iz) = xx.Values;
end
prob_vcam = npoint_vcam'/sum(npoint_vcam(:));


function [vcam, vcrm, nz, zint] = catVar(varname,timesteps)
vcam=[];vcrm=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D'},parm)
    zint = parm.zint/1000; % m to km
    nz = parm.nzi;
  end
  evalc(sprintf('vcam = cat(2,vcam,varCAM1D(1:parm.nzi,:,idv.varCAM1D.%s_cam));',varname));
  evalc(sprintf('vcrm = cat(2,vcrm,varCAM1D(1:parm.nzi,:,idv.varCAM1D.%s_crm));',varname));
end
