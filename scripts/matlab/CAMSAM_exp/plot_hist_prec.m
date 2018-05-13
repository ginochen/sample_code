% for plots under ~/archive/matlab/figure/precip_rate/
% scatterhist_allsizemcs_precc_precl_cam_spcam_tstep41-10-301.fig 
% scatterhist_precc_precl_mcc_cam_spcam_tstep41-10-401.fig
figure
timesteps = 41:10:401;
maxx=7;
vcam_conv=[];vcam_strat=[];
vcrm_conv=[];vcrm_strat=[];
for it = timesteps
  load(['var_PC1_' num2str(it) '.mat'],'varCRM0D','varCAM0D','parm')
  if it==timesteps(1)
    idv = varindex({'varCAM0D','varCRM0D'},parm);
  end
  vcam_conv = cat(1,vcam_conv,varCAM0D(:,idv.varCAM0D.precc));
  vcrm_conv = cat(1,vcrm_conv,varCRM0D(:,idv.varCRM0D.precc));
  vcam_strat = cat(1,vcam_strat,varCAM0D(:,idv.varCAM0D.precl));
  vcrm_strat = cat(1,vcrm_strat,varCRM0D(:,idv.varCRM0D.precl));
end
vcrm = [vcrm_conv; vcrm_strat];
vcam = [vcam_conv; vcam_strat];
vsize = size(vcrm_conv,1);
species = [1*ones(vsize,1); 2*ones(vsize,1)]; % convective (1), stratiform (2)
scatterhist(vcrm,vcam,'Group',species,'Color','rb','Marker','o+'); hold on
plot([0 maxx], [0 maxx],'k')
h = legend('Convective','Stratiform')
set(h,'fontsize',20,'box','off')
ylim([0 maxx])
xlim([0 maxx])
xlabel('P_{SPCAM} [mm/hr]','Fontsize',20)
ylabel('P_{CAM} [mm/hr]','Fontsize',20)
set(gca,'Fontsize',14,'box','off')
axis square
set(gcf,'color','w')
