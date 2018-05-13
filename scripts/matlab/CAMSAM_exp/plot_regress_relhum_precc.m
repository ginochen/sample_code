%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose : scatter plot the PRECT (surface precipitation rate) vs RELHUM for CAM, and fit a
% regress curve which must match up with the daily mean curve p =
% exp(15.6(r - 0.603)); and monthly mean curve p = exp(11.4(r -
% 0.522));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% do conditions
%%%%%%%%%%%%%%%%%%
docase = 2; % 1: CAM, 2: SPCAM
%%%%%%%%%%%%%%%%%%
% Param
%%%%%%%%%%%%%%%%%%
load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/index/basin_index.mat'); % iiocn, iinocn, basin
iilev= [14:30]; % selected vertical levels
rr=[0.2:0.01:0.85]; % relhum linear values for the regression curve
%%%%%%%%%%%%%%%%%%
% 3D var
%%%%%%%%%%%%%%%%%%
if (exist('idv','var')==0 | isfield(idv,'rh')==0 | isfield(idv,'rhcrm')==0) % if idv nonexist or rh field nonexist
   switch docase
   case 1
      idv.rh=3;     varname{idv.rh}    = 'RELHUM';   var{idv.rh}    = loadvar(varname{idv.rh},camCase,0.01,1);
   case 2
      idv.rhcrm=47; varname{idv.rhcrm} = 'RELHUM';   var{idv.rhcrm} = loadvar(varname{idv.rhcrm},spcamCase,0.01,0); % RELHUM for CRM
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (exist('idv','var')==0 | isfield(idv,'prect')==0 | isfield(idv,'spp_sfc'))
   switch docase    
   case 1
      idv.prect=5;    varname{idv.prect} = 'PRECT';    var{idv.prect} = loadvar(varname{idv.prect},camCase,spd*1000,1);% (mm/day) (convective + large scl) 
   case 2
      idv.spp_sfc=49; varname{idv.spp_sfc}   = 'SPPFLX_SFC';   var{idv.spp_sfc} = loadvar(varname{idv.spp_sfc},spcamCase,spd*1000,0); % Precip flux for CRM
%      idv.prectcrm=50;    varname{idv.prectcrm} = 'PRECT';    var{idv.prectcrm} = loadvar(varname{idv.prectcrm},camCase,spd*1000,1);% (mm/day) (convective + large scl) 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p-ave relhum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=1:numel(itimes)
   switch docase
   case 1
      relhum(:,:,it)=pwgtave(var{idv.rh}(:,:,:,itimes(it)),iilev,dlev);
   case 2
      relhum(:,:,it)=pwgtave(var{idv.rhcrm}(:,:,:,itimes(it)),iilev,dlev);
   end
end
relhumtmp= mean(relhum,3);
switch docase
case 1
   precctmp = mean(var{idv.prect},3); % time-ave
case 2
   precctmp = squeeze(var{idv.spp_sfc}(:,:,1,itimes));
%   precctmp = squeeze(var{idv.prectcrm}(:,:,itimes));
   precctmp(precctmp<0)=0;
   precctmp = mean(precctmp,3);
end
%
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear Regression Model 
%   ( for precipitation rate as a function of relative humidity 
%     log(precc) = b0 + b1*relhum
%     => precc = exp(b0+b1*relhum) )
for ib=1:numel(basin)+1
   precc20{ib}  = precctmp(iiocn{ib});
   relhum20{ib} = relhumtmp(iiocn{ib});
   logP20=log(precc20{ib});
   B{ib}=regress(logP20(:),[ones(numel(relhum20{ib}),1),relhum20{ib}(:)])
%   scatter(relhum20{ib}(:),precc20{ib}(:)); hold on; pause
end
figure
linetype = {'k:','k.','k--','k-o','r-'}; 
for ib=1:numel(basin)+1
   plot(rr,exp(B{ib}(1)+B{ib}(2).*rr),linetype{ib},'linewidth',2); hold on
end
legend([basin,'ALL'])
xlabel('Relative Humidity')
ylabel('P[mm/day]')
ylim([0 40])
switch docase
case 1
   title('CAM 14-day-mean P vs r')
  % save('prect_relhum_cam_20s20n_t700p14p30-ave.mat','B','ilat','ilon','iilev','relhum','iiocn','iinocn','basin','lat','lon','lev'); 
case 2
   title('SPCAM 14-day-mean P vs r')
  % save('sppflx_relhum_crm_20s20n_t700p14p30-ave.mat','B','ilat','ilon','iilev','relhum','iiocn','iinocn','basin','lat','lon','lev'); 
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% precc-Quantile-conditioned Regression
precc_quant = quantile(precc20(:),[0.25,0.5,0.75])
for i=1:numel(precc_quant)
   idx = find(precc20(:)<precc_quant(i));
   [B,BINT,R,RINT]=regress(logP20(idx),[ones(numel(idx),1),relhum20(idx)]);
   plot(rr,exp(B(1)+B(2).*rr),'k:','linewidth',2);     
   hold on
end
%[c, p, err, ind] = kmeans_clusters([precc20(:) relhum20(:)],3);
figure; hist(relhum20(:),50); xlim([0 1])
figure; hist(precc20(:),50)

