clear shdt zmdt
run ~/scripts/matlab/startup
main_stoch
idv.dt=1;   varname{idv.dt} = 'DTCOND'; var{idv.dt}   = loadvar(varname{idv.dt},  camCase,  spd,    0, [1,1,1,1], [Inf,Inf,Inf,500]);% (K/day)<---- only moist parametrized heating
idv.dcq=2;  varname{idv.dcq} = 'DCQ';   var{idv.dcq}  = loadvar(varname{idv.dcq}, camCase,  spd,    0, [1,1,1,1], [Inf,Inf,Inf,500]);% (kg/kg/day)<---- only moist parametrized water vapor tendency 
idv.spdq=3; varname{idv.spdq} = 'SPDQ'; var{idv.spdq} = loadvar(varname{idv.spdq},spcamCase,spd,    0, [1,1,1,1], [Inf,Inf,Inf,500]);% (kg/kg/day)<---- only moist parametrized water vapor tendency 
idv.spdt=30;varname{idv.spdt} = 'SPDT'; var{idv.spdt} = loadvar(varname{idv.spdt},spcamCase,spd,    0, [1,1,1,1], [Inf,Inf,Inf,500]);% T tendency due to CRM (without QRL+QRS radiative heating)
idv.pzm=10; varname{idv.pzm} = 'PRECCDZM';var{idv.pzm}= loadvar(varname{idv.pzm}, camCase,  spd*1e3,0, [1,1,1],   [Inf,Inf,500]);% (mm/day) 
idv.psh=11; varname{idv.psh} = 'PRECSH';  var{idv.psh}= loadvar(varname{idv.psh}, camCase,  spd*1e3,0, [1,1,1],   [Inf,Inf,500]);% (mm/day) 
idv.rh=3;   varname{idv.rh} = 'RELHUM';   var{idv.rh} = loadvar(varname{idv.rh},  camCase,  1e-2,   1, [1,1,1,1], [Inf,Inf,Inf,501]); % 501 equals 500 indices
idv.rhcrm=47;   varname{idv.rhcrm} = 'RELHUM';      var{idv.rhcrm}    = loadvar(varname{idv.rhcrm},   spcamCase, 0.01, 0, [1,1,1,1], [Inf,Inf,Inf,501] ); % RELHUM for CRM
idv.cmfmcdzm=6; varname{idv.cmfmcdzm} = 'CMFMCDZM'; var{idv.cmfmcdzm} = loadvar(varname{idv.cmfmcdzm},camCase,   1,    1, [1,1,1,1], [Inf,Inf,Inf,501]); % just read level from 850 to 985 skip 1000
idv.cmfmc=7;    varname{idv.cmfmc} = 'CMFMC';       var{idv.cmfmc}    = loadvar(varname{idv.cmfmc},   camCase,   1,    1, [1,1,1,1], [Inf,Inf,Inf,501]); % just read level from 850 to 985 skip 1000
% this code is used to find the basinwise deep/shallow convection heating
% profile, and decompose the modes
basin{1}='IND';  % indian ocean
cutoff_lat{1} = [[-20, 0],   [0,   20] ]; % lon=[40 to 120], lat=[-20 to 0], etc. cutoff_lat{1} = [-20, 0, 0, 20];
cutoff_lon{1} = [[40,  120], [40,  100]];
basin{2}='WPAC'; % western pacific
cutoff_lat{2} = [[-20, 0],   [0,   20] ];
cutoff_lon{2} = [[120, 180], [100, 180]];
basin{3}='EPAC'; % eastern pacific
cutoff_lat{3} = [[-20, 10],  [10,  20] ];
cutoff_lon{3} = [[180, 300], [180, 271]];
basin{4}='ATL';  % atlantic
cutoff_lat{4} = [[-20, 10],  [-20, 10],  [10,  20]];
cutoff_lon{4} = [[  0, 20],  [290, 360], [271, 350]];  
basin{5}='ALL';
%
landfrac = loadvar('LANDFRAC',camCase,1,0);
clear zmdt iillt_zm zmspdt
clear shdt iillt_sh shspdt 
for ib = 1:numel(basin)-1
tic
   iiizm=1; % counter for a basin-wise index (lon lat time) (shallow scheme)
   iiish=1; % counter for a basin-wise index (lon lat time) (zm scheme)
   for i = 1:numel(cutoff_lat{ib})/2 % divide by 2 gives the number of parts, since a basin is separated into two or three parts, and every part has 2 boundary index.
      iilon = find(lon>=cutoff_lon{ib}((i-1)*2+1) & lon<=cutoff_lon{ib}((i-1)*2+2))'; 
      iilat = find(lat>=cutoff_lat{ib}((i-1)*2+1) & lat<=cutoff_lat{ib}((i-1)*2+2))';
      for it=1:10:500
         for ii = iilon
            for jj = iilat
%               if (landfrac(ii,jj) == 0 & var{idv.pzm}(ii,jj,it)>0 & abs(var{idv.psh}(ii,jj,it))<=0.1 ) % if ocean index, if precip is deep
%               if (landfrac(ii,jj) == 0 & var{idv.pzm}(ii,jj,it)>0 ) % if ocean index, if precip is deep
%               if (landfrac(ii,jj) == 0 & var{idv.pzm}(ii,jj,it)>0 & var{idv.psh}(ii,jj,it)>0) % if ocean index, if precip is deep
                  %etmp{ib}(:,iii) = var{idv.ea}(ii,jj,:,it);
%                  zmdt{ib}(:,iiizm) = var{idv.dt}(ii,jj,:,it); % dQ for deep convection zm-scheme
%                  zmspdt{ib}(:,iiizm)= var{idv.spdt}(ii,jj,:,it); % dQ for SP conditioned on zm-scheme only (exclude shallow UW-scheme)
%                  zmdq{ib}(:,iiizm) = var{idv.dcq}(ii,jj,:,it); % dQ for deep convection zm-scheme
%                  zmspdq{ib}(:,iiizm)= var{idv.spdq}(ii,jj,:,it); % dQ for SP conditioned on zm-scheme only (exclude shallow UW-scheme)
%                  iillt_zm{ib}(:,iiizm) = [ii,jj,it]; % index of lon lat time when the condition is satisfied
%                  iiizm=iiizm+1;
%               end
               if (landfrac(ii,jj) == 0 & var{idv.psh}(ii,jj,it)>0 & abs(var{idv.pzm}(ii,jj,it))<=0.1) % if ocean index, if precip is deep
%                  %etmp{ib}(:,iii) = var{idv.ea}(ii,jj,:,it);
                  uwdq{ib}(:,iiish) = var{idv.dcq}(ii,jj,:,it); % dQ vertical profile at iiish index for deep convection zm-scheme
                  uwspdq{ib}(:,iiish)= var{idv.spdq}(ii,jj,:,it); % dQ for SP conditioned on zm-scheme only (exclude shallow UW-scheme)
%                  shdt{ib}(:,iiish) = var{idv.dt}(ii,jj,:,it); % dQ for shallow convection UW-scheme
%                  shspdt{ib}(:,iiish)= var{idv.spdt}(ii,jj,:,it); % dQ for SP conditioned on UW-scheme only (exclude (zm-scheme)
                  iillt_sh{ib}(:,iiish) = [ii,jj,it]; % index of lon lat time when the condition is satisfied
                  iiish=iiish+1;
               end
            end
         end
      end
   end
end
%
%
%
%
for ib = 1:numel(basin)-1
tic
   for i = 1:size(iillt_zm{ib},2)
      %rh_zm{ib}(:,i) = var{idv.rh}(iillt_zm{ib}(1,i),iillt_zm{ib}(2,i),:,iillt_zm{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
      cmfmc_zm{ib}(:,i) = var{idv.cmfmcdzm}(iillt_zm{ib}(1,i),iillt_zm{ib}(2,i),:,iillt_zm{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
      %rh_zmsp{ib}(:,i) = var{idv.rhcrm}(iillt_zm{ib}(1,i),iillt_zm{ib}(2,i),:,iillt_zm{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
   end
   for i = 1:size(iillt_sh{ib},2)
      cmfmc_sh{ib}(:,i) = var{idv.cmfmc}(iillt_sh{ib}(1,i),iillt_sh{ib}(2,i),:,iillt_sh{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
      %rh_sh{ib}(:,i) = var{idv.rh}(iillt_sh{ib}(1,i),iillt_sh{ib}(2,i),:,iillt_sh{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
      %rh_shsp{ib}(:,i) = var{idv.rhcrm}(iillt_sh{ib}(1,i),iillt_sh{ib}(2,i),:,iillt_sh{ib}(3,i)); % obtain all the deep convection indices from the saved index set iilt
   end
toc
end
%
%   
for ib = 1:numel(basin)-1
   uwddq{ib} = uwdq{ib} - uwspdq{ib};
   zmddq{ib} = zmdq{ib} - zmspdq{ib};
end
%
%
for ib = 1:numel(basin)-1
tic
%%%%% just use pca instead of svdVar !!!!!!!!!!!!!!!!!!!!!!!!!!!
   [svec_ddqzm{ib},sval_ddqzm{ib},pc_ddqzm{ib}] = svdVar(zmddq{ib});
%   [svec_dqzm{ib},sval_dqzm{ib},pc_dqzm{ib}] = svdVar(zmdq{ib});
%   [svec_dqzmsp{ib},sval_dqzmsp{ib},pc_dqzmsp{ib}] = svdVar(zmspdq{ib});
%   [svec_dquw{ib},sval_dquw{ib},pc_dquw{ib}] = svdVar(uwdq{ib});
%   [svec_dquwsp{ib},sval_dquwsp{ib},pc_dquwsp{ib}] = svdVar(uwspdq{ib});
   [svec_ddquw{ib},sval_ddquw{ib},pc_ddquw{ib}] = svdVar(uwddq{ib});
%   [svec_qzm{ib},sval_qzm{ib},pc_qzm{ib}] = svdVar(zmdt{ib});
%   [svec_qzmsp{ib},sval_qzmsp{ib},pc_qzmsp{ib}] = svdVar(zmspdt{ib});
%   [svec_cmz{ib},sval_cmz{ib},pc_cmz{ib}] = svdVar(cmfmc_zm{ib});
%   [svec_cms{ib},sval_cms{ib},pc_cms{ib}] = svdVar(cmfmc_sh{ib});
%   [svec_drhz{ib},sval_drhz{ib},pc_drhz{ib}] = svdVar(rh_zm{ib}-rh_zmsp{ib});
%   [svec_drhs{ib},sval_drhs{ib},pc_drhs{ib}] = svdVar(rh_sh{ib}-rh_shsp{ib});
%   [svec_rhz{ib},sval_rhz{ib},pc_rhz{ib}] = svdVar(rh_zm{ib});
%   [svec_rhs{ib},sval_rhs{ib},pc_rhs{ib}] = svdVar(rh_sh{ib});
%   [svec_dqz{ib},sval_dqz{ib},pc_dqz{ib}] = svdVar(zmspdt{ib}-zmdt{ib});
%   [svec_dqs{ib},sval_dqs{ib},pc_dqs{ib}] = svdVar(shspdt{ib}-shdt{ib});
%   [svec_z{ib},sval_z{ib},pc_z{ib}] = svdVar(zmdt{ib});
%   [svec_sps{ib},sval_sps{ib},pc_sps{ib}] = svdVar(shspdt{ib});
%   [svec_s{ib},sval_s{ib},pc_s{ib}] = svdVar(shdt{ib});
%   [svec_spz{ib},sval_spz{ib},pc_spz{ib}] = svdVar(zmspdt{ib});
   %[svec_sh{ib},sval_sh{ib},pc_sh{ib}] = svdVar(shdt{ib});
   %[svec_zm{ib},sval_zm{ib},pc_zm{ib}] = svdVar(zmdt{ib});
toc
end
% /projects/rsmas/kirtman/gchen/archive/matlab/figure/
%save('dtcond_spdt_ZMincUW_UWexcZM_svd-basinwise_t-1.mat','svec_z','sval_z','pc_z','svec_sps','sval_sps','pc_sps')         
%save('dtcond_spdt_ZMincUW_UWexcZM_svd-basinwise_t-1.mat','svec_s','sval_s','pc_s','svec_spz','sval_spz','pc_spz','-append')         
%save('dtcond_spdt_pzm_excludesh_psh_excludezm_svd-basinwise.mat','svec_s','sval_s','pc_s','svec_z','sval_z','pc_z')         
%save('dtcond_spdt_pzm_excludesh_psh_excludezm_svd-basinwise.mat','svec_sps','sval_sps','pc_sps','svec_spz','sval_spz','pc_spz','-append')
%save('dtcond_spdt_pzm_includesh_psh_excludezm_svd-basinwise.mat',...
%     'svec_spz','sval_spz','pc_spz','svec_z','sval_z','pc_z','lev','iillt_zm',...
%     'zmdt','zmspdt','cutoff_lat','cutoff_lon','basin')
%save('dQ_pzm_excludesh_psh_excludezm_svd-basinwise.mat',... <------ dtcond - spdt
%     'svec_dqz','sval_dqz','pc_dqz','lev',...
%     'svec_dqs','sval_dqs','pc_dqs',...
%     'cutoff_lat','cutoff_lon','basin')
%save('dRH_pzm_excludesh_psh_excludezm_svd-basinwise.mat',...
%     'svec_drhz','sval_drhz','pc_drhz','lev',...
%     'svec_drhs','sval_drhs','pc_drhs',...
%     'cutoff_lat','cutoff_lon','basin')
%save('RH_pzm_excludesh_psh_excludezm_svd-basinwise.mat',...
%     'svec_rhz','sval_rhz','pc_rhz','lev',...
%     'svec_rhs','sval_rhs','pc_rhs',...
%     'cutoff_lat','cutoff_lon','basin')
%save('RH_pzm_excludesh_psh_excludezm_svd-basinwise.mat','rh_zm','rh_sh','-append')
%save('CMFMC_pzm_excludesh_psh_excludezm_svd-basinwise.mat',...
%     'svec_cmz','sval_cmz','pc_cmz','lev',...
%     'svec_cms','sval_cms','pc_cms',...
%     'cutoff_lat','cutoff_lon','basin')
%save('CMFMC_pzm_excludesh_psh_excludezm_svd-basinwise.mat','cmfmc_zm','cmfmc_sh','-append')
%save('dtcond_spdt_ZMintersecUW_svd-basinwise.mat',...
%     'svec_qzm','sval_qzm','pc_qzm','lev',...
%     'svec_qzmsp','sval_qzmsp','pc_qzmsp',...
%     'cutoff_lat','cutoff_lon','basin')
%save('dq_spdq_ZMexcUW_svd-basinwise.mat',...
%     'svec_dqzm','sval_dqzm','pc_dqzm','lev',...
%     'svec_dqzmsp','sval_dqzmsp','pc_dqzmsp',...
%     'cutoff_lat','cutoff_lon','basin')
%save('dq_spdq_UWexcZM_svd-basinwise.mat',...
%     'svec_dquw','sval_dquw','pc_dquw','lev',...
%     'svec_dquwsp','sval_dquwsp','pc_dquwsp',...
%     'cutoff_lat','cutoff_lon','basin')
save('ddq_ZMexcUW_UWexcZM_svd-basinwise.mat',...
     'svec_ddquw','sval_ddquw','pc_ddquw','lev',...
     'svec_ddqzm','sval_ddqzm','pc_ddqzm',...
     'cutoff_lat','cutoff_lon','basin')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot basin-wise mean Q projected on a EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
im=4; % im-th EOF 
for ib=1:numel(basin)-1
   subplot(1,4,ib)
   hold on; 
   plot(lev,svec_z{ib}(:,im)*mean(pc_z{ib}(:,im))','b');  % Q_CAM projected on EOF-im
   plot(lev,svec_spz{ib}(:,im)*mean(pc_spz{ib}(:,im))','r'); % Q_SPCAM projected on EOF-im
   view(90,90); % rotate the plot to vertical pressure level
   ylim([-6 6]);
   title(basin{ib}) % give basin name to subplot title
end
sh = subplot(1,4,4) % define the handle for subplot to use for legend
lh = legend(sh,'Q_{CAM}','Q_{SPCAM}') % lengend on the last panel
set(lh,'position',[0.9 0.9 .1 .1]) % set the legend position on the top-right corner
suptitle(['Time mean Q projected on EOF' num2str(im)]) % give a supertitle over all subplots
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot tropical ocean (basin-wise average) mean Q projected on a EOF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for im=1
%   subplot(1,2,1)
   hold on
   cat_pc_rhz = [];
   for ib = 1:4; cat_pc_rhz = horzcat(cat_pc_rhz,svec_rhz{ib}(:,im)*pc_rhz{ib}(:,im)'); end % concat all basin into one 30 lev * # of space-time matrix
   rh_eof_std_zm = std(cat_pc_rhz,0,2);
   rh_eof_ave_zm = mean(cat_pc_rhz,2);
   plot(lev,rh_eof_ave_zm+rh_eof_std_zm,'r:','linewidth',1.5);
   plot(lev,rh_eof_ave_zm-rh_eof_std_zm,'r:','linewidth',1.5);
   plot(lev,rh_eof_ave_zm,'r','linewidth',1.5);
%   view(90,90)
%   subplot(1,2,2)
   hold on
   cat_pc_rhs = [];
   for ib = 1:4; cat_pc_rhs = horzcat(cat_pc_rhs,svec_rhs{ib}(:,im)*pc_rhs{ib}(:,im)'); end % concat all basin into one 30 lev * # of space-time matrix
   rh_eof_std_sh = std(cat_pc_rhs,0,2);
   rh_eof_ave_sh = mean(cat_pc_rhs,2);
   plot(lev,rh_eof_ave_sh+rh_eof_std_sh,'k:','linewidth',1.5);
   plot(lev,rh_eof_ave_sh-rh_eof_std_sh,'k:','linewidth',1.5);
   plot(lev,rh_eof_ave_sh,'k','linewidth',1.5);
   view(90,90)
end
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot tropical ocean (basin-wise average) mean Q projected on a EOF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure;
set(h,'position',[22         196        1376         606]);
nb = numel(basin)-1; % # of basins
icase = 9; % 1) ZM 2) UW 3) dQZM & dQUW 
nm=4; %total numer of EOF modes
for im=1:nm % EOF modes
   subplot(1,nm,im)
   cat_pc_spdtz=[]; cat_pc_dtz=[]; cat_pc_spdts=[]; cat_pc_dts=[]; cat_pc_dQz=[]; cat_pc_dQs=[]; cat_pc_cmz=[]; cat_pc_cms=[];
   cat_pc_qzm=[]; cat_pc_qzmsp=[];
   cat_pc_dqzmsp=[]; cat_pc_dqzm=[]; cat_pc_ddqzm=[]; cat_pc_dquwsp=[]; cat_pc_dquw=[]; cat_pc_ddquw=[];
   for ib=1:numel(basin)-1
      if icase == 1
         cat_pc_spdtz = horzcat(cat_pc_spdtz, svec_spz{ib}(:,im)*pc_spz{ib}(:,im)');
         cat_pc_dtz = horzcat(cat_pc_dtz,svec_z{ib}(:,im)*pc_z{ib}(:,im)');
      elseif icase == 2
         cat_pc_spdts = horzcat(cat_pc_spdts,svec_sps{ib}(:,im)*pc_sps{ib}(:,im)');
         cat_pc_dts = horzcat(cat_pc_dts,svec_s{ib}(:,im)*pc_s{ib}(:,im)');
      elseif icase == 3
         cat_pc_dQz = horzcat(cat_pc_dQz,svec_dqz{ib}(:,im)*pc_dqz{ib}(:,im)');
         cat_pc_dQs = horzcat(cat_pc_dQs,svec_dqs{ib}(:,im)*pc_dqs{ib}(:,im)');
      elseif icase == 4
         cat_pc_cmz = horzcat(cat_pc_cmz,svec_cmz{ib}(:,im)*pc_cmz{ib}(:,im)');
         cat_pc_cms = horzcat(cat_pc_cms,svec_cms{ib}(:,im)*pc_cms{ib}(:,im)');
      elseif icase == 5
         cat_pc_qzm   = horzcat(cat_pc_qzm,svec_qzm{ib}(:,im)*pc_qzm{ib}(:,im)');
         cat_pc_qzmsp = horzcat(cat_pc_qzmsp,svec_qzmsp{ib}(:,im)*pc_qzmsp{ib}(:,im)');
      elseif icase == 6
         cat_pc_spdtz = horzcat(cat_pc_spdtz, svec_spz{ib}(:,im)*pc_spz{ib}(:,im)');
         cat_pc_dtz = horzcat(cat_pc_dtz,svec_z{ib}(:,im)*pc_z{ib}(:,im)');
         cat_pc_dQz = horzcat(cat_pc_dQz,svec_dqz{ib}(:,im)*pc_dqz{ib}(:,im)');
      elseif icase == 7
         cat_pc_spdts = horzcat(cat_pc_spdts,svec_sps{ib}(:,im)*pc_sps{ib}(:,im)');
         cat_pc_dts = horzcat(cat_pc_dts,svec_s{ib}(:,im)*pc_s{ib}(:,im)');
         cat_pc_dQs = horzcat(cat_pc_dQs,svec_dqs{ib}(:,im)*pc_dqs{ib}(:,im)');
      elseif icase == 8
         cat_pc_dqzmsp = horzcat(cat_pc_dqzmsp,svec_dqzmsp{ib}(:,im)*pc_dqzmsp{ib}(:,im)');
         cat_pc_dqzm = horzcat(cat_pc_dqzm,svec_dqzm{ib}(:,im)*pc_dqzm{ib}(:,im)');
         cat_pc_ddqzm = horzcat(cat_pc_ddqzm,svec_ddqzm{ib}(:,im)*pc_ddqzm{ib}(:,im)');
      elseif icase == 9
         cat_pc_dquwsp = horzcat(cat_pc_dquwsp,svec_dquwsp{ib}(:,im)*pc_dquwsp{ib}(:,im)');
         cat_pc_dquw   = horzcat(cat_pc_dquw,svec_dquw{ib}(:,im)*pc_dquw{ib}(:,im)');
         cat_pc_ddquw  = horzcat(cat_pc_ddquw,svec_ddquw{ib}(:,im)*pc_ddquw{ib}(:,im)');
      else
         error('specify to do ZM-deep or UW-shallow condition');
      end
   end
   hold on;
   if (icase == 1 )
      plot(lev,mean(cat_pc_dtz,2),'r','linewidth',1.5);
      plot(lev,mean(cat_pc_spdtz,2),'color',[.5 .5 .5],'linewidth',1.2);
   elseif (icase == 2)
      plot(lev,mean(cat_pc_dts,2),'r','linewidth',1.5);
      plot(lev,mean(cat_pc_spdts,2),'color',[.5 .5 .5],'linewidth',1.2);
   elseif (icase == 3)
      plot(lev,-mean(cat_pc_dQz,2),'r','linewidth',2); %negative sign is for dQ to become Q_CAM - Q_SPCAM
      plot(lev,-mean(cat_pc_dQs,2),'color',[.5 .5 .5],'linewidth',2);
   elseif (icase == 4)
      plot(ilev,mean(cat_pc_cmz,2),'r','linewidth',1.5); %negative sign is for dQ to become Q_CAM - Q_SPCAM
      plot(ilev,mean(cat_pc_cms,2),'color',[.5 .5 .5],'linewidth',1.2);
%      plot(lev,RHeof_basin_sumzm/nb,'r','linewidth',1.5); 
%      plot(lev,RHeof_basin_sumsh/nb,'color',[.5 .5 .5],'linewidth',1.2);
%      for it = 1:50:size(pc_rhz{ib},1); plot(lev,svec_rhz{ib}(:,im)*pc_rhz{ib}(it,im),'color',[.5 .5 .5]);hold on; end
   elseif (icase == 5)
      plot(lev,mean(cat_pc_qzm,2),'r','linewidth',1.5); %negative sign is for dQ to become Q_CAM - Q_SPCAM
      plot(lev,mean(cat_pc_qzmsp,2),'color',[.5 .5 .5],'linewidth',1.2);
   elseif icase == 6
      plot(lev,mean(cat_pc_dtz,2),'k','linewidth',2);
      plot(lev,mean(cat_pc_spdtz,2),'color',[.5 .5 .5],'linewidth',2);
      plot(lev,-mean(cat_pc_dQz,2),'r:','linewidth',2); %negative sign is for dQ to become Q_CAM - Q_SPCAM
   elseif icase == 7
      plot(lev,mean(cat_pc_dts,2),'k','linewidth',2);
      plot(lev,mean(cat_pc_spdts,2),'color',[.5 .5 .5],'linewidth',2,'linestyle','--');
      plot(lev,-mean(cat_pc_dQs,2),'r:','linewidth',2);
   elseif (icase == 8 )
      plot(lev,mean(cat_pc_dqzm,2),'k','linewidth',2);
      plot(lev,mean(cat_pc_dqzmsp,2),'color',[.5 .5 .5],'linewidth',2,'linestyle','--');
      plot(lev,-mean(cat_pc_ddqzm,2),'r:','linewidth',2);
   elseif (icase == 9 )
      plot(lev,mean(cat_pc_dquw,2),'k','linewidth',2);
      plot(lev,mean(cat_pc_dquwsp,2),'color',[.5 .5 .5],'linewidth',2,'linestyle','--');
      plot(lev,-mean(cat_pc_ddquw,2),'r:','linewidth',2);
   end
   view(90,90);
   if im == 1; xlabel('Pressure (hPa)','Fontsize',16); end
   if icase == 1 | icase == 2 | icase == 5 | icase == 6 | icase == 7
      ylim([-5.1 5.1]);
   elseif icase == 8 | icase == 9
      ylim([-0.03 0.03])
   elseif icase == 3
      ylim([-6 6])
   elseif icase == 4
      if im == 1; ylim([-0.005 0.05]); else; ylim([-0.005 0.005]);    end
   end
   title(['EOF' num2str(im)],'Fontsize',20)
end
%
%
sh = subplot(1,4,nm); % define the handle for subplot to use for legend
if icase ==1 | icase ==2
   lh = legend(sh,'Q_{CAM}','Q_{SPCAM}'); % lengend on the last panel
elseif icase == 3
   lh = legend(sh,'dQ_{deep}','dQ_{shallow}'); % lengend on the last panel
elseif icase == 4
   lh = legend(sh,'M_{deep}','M_{shallow}'); % lengend on the last panel
elseif icase == 6 | icase ==7
   lh = legend(sh,'Q_{CAM}','Q_{SPCAM}','dQ');
end
set(lh,'position',[0.88 0.87 .1 .1]); % set the legend position on the top-right corner
set(lh,'Fontsize',16);
set(lh,'EdgeColor',[1 1 1])
if icase == 1;
   sth = suptitle(['Tropical ocean time-mean Deep Q projected on EOFs']); % give a supertitle over all subplots
elseif icase == 2;
   sth = suptitle(['Tropical ocean time-mean Shallow Q projected on EOFs']); % give a supertitle over all subplots
elseif icase == 3;
   sth =  suptitle(['Tropical ocean time-mean dQ=Q_{CAM}-Q_{SPCAM} projected on EOFs']); % give a supertitle over all subplots
elseif icase == 4;
   sth =  suptitle(['Tropical ocean time-mean M_u projected on EOFs']); % give a supertitle over all subplots
elseif icase == 5;
   sth =  suptitle(['Tropical ocean time-mean Deep Q projected on EOFs']); % give a supertitle over all subplots
elseif icase == 6;
else
   error('specify to do ZM-deep or UW-shallow condition');
end
set(sth,'Fontsize',18);
if icase <=2;
   slh = suplabel('Q (K/day)','x',[0.45 0.12 0.1 0.1]);
elseif icase == 3
   slh = suplabel('dQ (K/day)','x',[0.45 0.12 0.1 0.1]);
elseif icase == 4
   slh = suplabel('Mass Flux (kg/s m^{-2})','x',[0.45 0.1 0.1 0.1]);
end
set(slh,'Fontsize',16);
set(h,'color',[1 1 1]);

% ilat = find(lat>=-15 & lat<=15);
% dtcond = var{idv.dt}(:,ilat,:,1);
% [ix,iy] = find(var{idv.psh}(:,ilat,1)>0);
% ii=1; for i=1:numel(ix); shdt(ii,:) = dtcond(ix(i),iy(i),:); ii=ii+1; end
% [svec_sh,sval_sh,pc_sh] = svdVar(shdt');

 % plot sh heating top 3 modes (time-ave)
 figure; hold on
 j=1;plot(lev,mean(pc_sh{ib}(:,j))*sval_sh{ib}((j)*svec_sh{ib}((:,j),'r','linewidth',2)
 j=2;plot(lev,mean(pc_sh{ib}((:,j))*sval_sh{ib}((j)*svec_sh{ib}((:,j),'k','linewidth',2)  
 j=3;plot(lev,mean(pc_sh{ib}((:,j))*sval_sh{ib}((j)*svec_sh{ib}((:,j),'k--','linewidth',2)

% [ixz,iyz] = find(var{idv.pzm}(:,ilat,1)>0);
% ii=1; for i=1:numel(ixz); zmdt(ii,:) = dtcond(ixz(i),iyz(i),:); ii=ii+1; end
% [svec_zm,sval_sh,pc_zm] = svdVar(zmdt');

 % plot zm heating top 3 modes (time-ave)
 figure; hold on
 j=1;plot(lev,mean(pc_zm{ib}((:,j))*sval_zm{ib}((j)*svec_zm{ib}((:,j),'r','linewidth',2)
 j=2;plot(lev,mean(pc_zm{ib}((:,j))*sval_zm{ib}((j)*svec_zm{ib}((:,j),'k','linewidth',2)  
 j=3;plot(lev,mean(pc_zm{ib}((:,j))*sval_zm{ib}((j)*svec_zm{ib}((:,j),'k--','linewidth',2)
