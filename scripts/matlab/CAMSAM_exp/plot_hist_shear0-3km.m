% for plots under ~/archive/matlab/figure/shear/
% scatterhist_allsizemcs_shear_cam_spcam_3km_10km_tstep41-10-301.fig
% scatterhist_mcc_shear_cam_spcam_3km_10km_tstep41-10-401.fig
maxx=35;
vcam_du3km = []; vcrm_du3km = []; vcam_du10km = []; vcrm_du10km = [];
timesteps=41:10:401;
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'stabilityIndex','parm')
  vcam_du3km = cat(1,vcam_du3km,stabilityIndex(:,3));
  vcrm_du3km = cat(1,vcrm_du3km,stabilityIndex(:,4));
  vcam_du10km = cat(1,vcam_du10km,stabilityIndex(:,5));
  vcrm_du10km = cat(1,vcrm_du10km,stabilityIndex(:,6));
end
vcam = [vcam_du3km*parm.dz3km; vcam_du10km*parm.dz10km]; 
vcrm = [vcrm_du3km*parm.dz3km; vcrm_du10km*parm.dz10km];
vsize = size(vcam_du3km,1);
species = [1*ones(vsize,1); 2*ones(vsize,1)];
scatterhist(vcrm,vcam,'Group',species,'Color','br','Marker','o.'); hold on
plot([0 maxx], [0 maxx],'k')
h = legend('|U_{3km}-U_{0km}|','|U_{10km}-U_{0km}|')
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
