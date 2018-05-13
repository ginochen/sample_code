% for plots under ~/archive/matlab/figure/shear/
% scatterhist_allsizemcs_shear_cam_spcam_3km_10km_tstep41-10-301.fig
% scatterhist_mcc_shear_cam_spcam_3km_10km_tstep41-10-401.fig
maxx=55; % max shear for both x and y axis
timesteps=11:10:401;
load('var_PC1_MCSt.mat','varCAM0D_mcsall','varCRM0D_mcsall')
load('var_PC1_41.mat','parm'); % load parm to get the variable indices assoc to varCAM0D

idv = getVarIndex({'varCAM0D','varCRM0D'},parm);
vcam_du200850 = varCAM0D_mcsall(idv.varCAM0D.du200850mb,:)';
vcam_du03km = varCAM0D_mcsall(idv.varCAM0D.du03km,:)';
vcam_du06km = varCAM0D_mcsall(idv.varCAM0D.du06km,:)';
vcam_du010km = varCAM0D_mcsall(idv.varCAM0D.du010km,:)';
vcrm_du200850 = varCRM0D_mcsall(idv.varCRM0D.du200850mb,:)';
vcrm_du03km = varCRM0D_mcsall(idv.varCRM0D.du03km,:)';
vcrm_du06km = varCRM0D_mcsall(idv.varCRM0D.du06km,:)';
vcrm_du010km = varCRM0D_mcsall(idv.varCRM0D.du010km,:)';

vsize = size(vcam_du200850,1);
vcam = [vcam_du03km; vcam_du06km; vcam_du010km];
vcrm = [vcrm_du03km; vcrm_du06km; vcrm_du010km];
species = [1*ones(vsize,1); 2*ones(vsize,1); 3*ones(vsize,1)];
scatterhist(vcrm,vcam,'Group',species,'Color','rbk','Marker','o.^'); hold on
plot([0 maxx], [0 maxx],'k')
h = legend('|U_{0.25km}-U_{3km}|','|U_{0.25km}-U_{6km}|','|U_{0.25km}-U_{10km}|')
set(h,'fontsize',20,'box','off')
ylim([0 maxx])
xlim([0 maxx])
xlabel('dU_{SPCAM} [m/s]','Fontsize',20)
ylabel('dU_{CAM} [m/s]','Fontsize',20)
set(gca,'Fontsize',14,'box','off')
axis square
set(gcf,'color','w')
xticks([ 10,20,30]);
yticks([ 10,20,30]);
%save('~/archive/matlab/figure/shear/scatterhist_mcsClusters_g15hrs_shear_cam_spcam_3km_10km_tstep11-10-301.fig')
