function ha = plot_waterphase(it,iillt,varCAM0D,varCAM1D,varCRM0D,varCRM1D,varCRM2D,parm,stabilityIndex,mcs,arxv)
%function ha = plot_waterphase(it,iillt,varCAM0D,varCAM1D,varCRM0D,varCRM1D,varCRM2D,parm,stabilityIndex,mcs,arxv)
% var and varCAM are both interpolated onto height coordinate (200m to 18km)
% plot basinwise spatially averaged organized cloud water over lead time
% plot TS,LFLX in rh0 to show that the difference of strength is due
% to SST
% plot EOF1 vert profile for each basin at the side
%
% Conclusions: 1) spcam has low level lcwat peak, which may moisten
% the low level and favor future convection (Tompkins 2001), 
% 2) cam only has upper level lcwat peak in most cases, so a drier low
% level atmosphere
% 
% Plots:
% Extreme cases (strong deep conv & organization conditioned on MCS
% index)
% 1) Mixing ratio time lapse for extreme cases + rainrates (condi
% 2) organization time-scale 
% Comparison of different cloud types in SPCAM vs CAM
% 
% parm.lagind*30min
%
display(mcs.units(iillt))
plot_crmQ = 1;
plot_LSQ = 0;
plot_heating = 0;
plot_crmrain = 1;
plot_relhum = 0;
plot_stability = 0;
str = {'varCAM0D','varCAM1D','varCRM0D','varCRM1D','varCRM2D','stabilityIndex'};
idv = getVarIndex(str,parm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% turn CIN NaN values to zero %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stabilityIndex(isnan(stabilityIndex))=0;
run ~/scripts/matlab/startup
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  determine the archive path and load variable %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if ( isempty(var) || isempty(varCAM) || isempty(parm) || isempty(stabilityIndex))
%   if (arxv) % 1 for working on macbook 0 for working on pegasus
%      archivePath = '/Users/g/archive/matlab/var/organizedQcQiQr/18km_lons11_lead0-10-lead48_actual/'
%   else % if working on pegasus
%      archivePath = '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/18km_lons11_lead0-10-lead48_actual/';
%   end
%   load([archivePath 'var_PC1_' num2str(it) '.mat'],'var','varCAM','stabilityIndex','parm')
%end
%  
%
condm = input('mean (1), median (2)');
condkm  = input('kmean (1), kmedoid (2), enter if do nothing'); 
if ( exist(condkm)==0 ); 
   condkm=0; 
elseif ( ismember(condkm,[1,2])  ) 
   kc = 2; % number of clusters, one cluster kmean is just the regular mean
   maxiter = 1000; % maximum number of interations before convergence for kmeans
   if ( condkm == 1 ) 
      kmean_medoid = 'kmeans'; 
   elseif ( condkm == 2 ) 
      kmean_medoid = 'kmedoids'; 
   end
end % 'median'
% mean shows double peak (due to averaging outliers) for lcwat in upper and lower atm, wereas median shows the realistic peak at top or bottom only
if condm == 1; mean_median = 'mean'; nanmean_nanmedian = 'nanmean'; elseif condm == 2; mean_median = 'median'; nanmean_nanmedian='nanmedian'; end % 'median'
qthres=[0.01 0.01]; 
%nlons=size(var{1},1); 
nlons=parm.nlons;
nt=parm.nlagt; %size(var{1},4); 
%
zint = parm.zint/1000; % in km
[X Z]=meshgrid([1:parm.nx],zint);
for ib=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lon-lat indices conditioned on stability %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   if isempty(iillt)
%      ii11=[];ii22=[];
%      for iilon = 1:nlons
%         % find deep convection before shear kicks in 
%         ii1 = find(squeeze(stabilityIndex(iilon,1,:,idv.LI_cam))<0.051); % LI is from MCS index
%         ii11 = cat(1,ii1,ii11);
%         ii2 = find((squeeze(stabilityIndex(iilon,end,:,idv.du3km_cam)))*1000>3); %0-3km shear is from MCS index
%         ii22 = cat(1,ii2,ii22);
%      end
%      %ii2 = find(squeeze(stabilityIndex(6,5,:,idv.cape_cam))>1000)
%      iillt = intersect(ii11,ii22)
%      iillt = ii11; 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % convective or stratiform %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %cloudtype(iilon,:,itlag,
%      if (isempty(iillt))
%         break
%      end
%   end
   tic 
   nsamp = size(varCRM2D,5); % number of big PC1 spatial samples
   figure(2*(ib-1)+1); 
   [ha ~] = tight_subplot(nt,nlons,[0.02,0],[.05 .05],[.05 .05]); %break
   iii=1; % subplot counter
   for iit=1:nt; % the lag time indices
      %%%%%%%%%%%%%%%%%%%%%%%%%
      % Define Variables here %
      %%%%%%%%%%%%%%%%%%%%%%%%%
      % kmeans or kmedoids
      if ( ismember(condkm,[1,2]) ) % kmeans or kmedoids
         qc_crm2D_km     = kmeans_kmedoids3D(idv.qc,    1000, var,    ib, iit, iillt, nsamp, parm, kc, kmean_medoid);
         qi_crm2D_km     = kmeans_kmedoids3D(idv.qi,    1000, var,    ib, iit, iillt, nsamp, parm, kc, kmean_medoid);
         qr_crm2D_km     = kmeans_kmedoids3D(idv.qr,    1000, var,    ib, iit, iillt, nsamp, parm, kc, kmean_medoid);
         lcwat_cam1D_km  = kmeans_kmedoids2D(idv.lcwat, 1000, varCAM, ib, iit, iillt, nsamp, parm, kc, kmean_medoid, maxiter);
         dtcond_cam1D_km = kmeans_kmedoids2D(idv.dtcond,   1, varCAM, ib, iit, iillt, nsamp, parm, kc, kmean_medoid, maxiter);
         spdt_cam1D_km   = kmeans_kmedoids2D(idv.spdt,     1, varCAM, ib, iit, iillt, nsamp, parm, kc, kmean_medoid, maxiter);
      end
      % mean or median
      qc_crm2D_m      = mean_median2D(idv.varCRM2D.qc, 1000, varCRM2D,    ib, iit, iillt, mean_median); 
      qi_crm2D_m      = mean_median2D(idv.varCRM2D.qi, 1000, varCRM2D,    ib, iit, iillt, mean_median); 
      lcwat_crm2D_m   = qc_crm2D_m + qi_crm2D_m; % this is how CAM defines, cldice + cldliq 
      qr_crm2D_m      = mean_median2D(idv.varCRM2D.qr, 1000, varCRM2D,    ib, iit, iillt, mean_median);
      rh_crm2D_m      = mean_median2D(idv.varCRM2D.rh,    1, varCRM2D,    ib, iit, iillt, mean_median); 
      u_crm2D_m       = mean_median2D(idv.varCRM2D.u,     1, varCRM2D,    ib, iit, iillt, mean_median);
      w_crm2D_m       = mean_median2D(idv.varCRM2D.w,     1, varCRM2D,    ib, iit, iillt, mean_median);
      spdt_cam1D_m    = mean_median1D(idv.varCAM1D.spdt, 1,    varCAM1D, ib, iit, iillt, mean_median); % [K/day]
      dtcond_cam1D_m  = mean_median1D(idv.varCAM1D.dtcond, 1,    varCAM1D, ib, iit, iillt, mean_median); % [K/day]
      lcwat_cam1D_m   = mean_median1D(idv.varCAM1D.lcwat, 1000, varCAM1D, ib, iit, iillt, mean_median); % [g/kg]
      relhum_cam1D_m  = mean_median1D(idv.varCAM1D.relhum, 1, varCAM1D, ib, iit, iillt, mean_median); % [g/kg]
      lcwat_crm1D_m   = mean_median1D([],[], lcwat_crm2D_m, [], [], [], mean_median); 
      rh_crm1D_m      = mean_median1D([],[], rh_crm2D_m,    [], [], [], mean_median); 
%      cape            = mean_median0D([idv.cape_cam,idv.cape_crm], 1, stabilityIndex, ib, iit, iillt, mean_median); 
%      cin             = mean_median0D([idv.cin_cam,idv.cin_crm], 1, stabilityIndex, ib, iit, iillt, mean_median); 
      LI              = mean_median0D([idv.stabilityIndex.LI_cam,idv.stabilityIndex.LI_crm], 1, stabilityIndex, ib, iit, iillt, mean_median); 
      cloudtype       = mean_median1D(idv.varCRM1D.cloudtype, 1, varCRM1D, ib, iit, iillt, nanmean_nanmedian);
      rainrate        = mean_median1D(idv.varCRM1D.rainrate,  1, varCRM1D, ib, iit, iillt, mean_median);
%      diffdt_m        = diffmean_median1D(3,2, 1, varCAM1D, ib, iit, iillt, mean_median);
      diffdt_m        = spdt_cam1D_m  - dtcond_cam1D_m; % [K/day]
      difflcwat_m     = lcwat_crm1D_m - lcwat_cam1D_m; % [g/kg]
      %%%%%%%%%%%%%%%%%%%%%
      %  Heating Contour  %
      %%%%%%%%%%%%%%%%%%%%%
%      figure(2*(ib-1)+1); subplot(nt,1,iit)
%      contourf(1:nlons,zint,diffdt_m',50,'linestyle','none'); caxis([0 50])
%      figure(2*(ib-1)+2); subplot(nt,1,iit)
%      contourf(1:nlons,zint,difflcwat_m',50,'linestyle','none'); caxis([0 0.1])
      for iil = 1:nlons
%break
         figure(2*(ib-1)+1);
         %%%%%%%%%%%%%%%%%%%
         % Start plot here %
         %%%%%%%%%%%%%%%%%%%
         ax = ha(iii); ax1p = ax.Position; 
         %  
         %%%%%%%%%%%%%%%%%%%%%%%
         %    QICE+QLIQ+QRAIN  %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_crmQ)
         ax(1) = axes('Position',ax1p,'Ylim',[zint(1) zint(end)]);% AXES NEEDS TO BE CALLED BEFORE PLOTTING!!! 
         %h1 = contourf(X, Z, lcwat_crm(iil,:,:)', qthres, 'linestyle', 'none', 'facecolor', [.7 .7 .7], 'Parent',ax(1)); hold on
         h1 = contourf(X, Z, squeeze(lcwat_crm2D_m(iil,:,:,1))', qthres, 'linestyle', 'none', 'facecolor', [.7 .7 .7], 'Parent',ax(1)); hold on
         h2 = contourf(X, Z, squeeze(qr_crm2D_m(iil,:,:,1))',  qthres, 'linestyle', 'none', 'facecolor', [.5 .5 .5], 'Parent',ax(1)); hold on
         %h2 = contourf(X, Z, qr_crm,  qthres, 'linestyle', 'none', 'facecolor', [.5 .5 .5], 'Parent',ax(1)); 
         h3 = streamslice(X,Z,squeeze(u_crm2D_m(iil,:,:))',squeeze(w_crm2D_m(iil,:,:))','Parent',ax(1));  set(h3,'color',[.6 .6 .6]); 
         end
         % plot lcwat
         %
         %%%%%%%%%%%%%%%%%%%%%%%
         %   LS QICE+QLIQ      %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_LSQ)
         if iii <= 11;  xtickval2 = [ -.03 0 .03]; xtickval3 = [-10.5 0 10.5]; xtickval4=[0 20]; else; xtickval2 = ''; xtickval3=''; end
         ax(2) = axes('Position',ax1p,'YTickLabel','','XAxisLocation','top',...
                      'Xcolor','r','Color','none','Xlim',[-0.05 0.05],'Ylim',[zint(1) zint(end)],...
                      'XTick',xtickval2); hold on % <------- Make sure this hold on is here, otherwise it'll be overwritten
%         h4 = plot(difflcwat_m(iil,:,1),       zint, 'r-', 'Parent', ax(2)); 
         h4 = plot(lcwat_cam1D_m(iil,:,1),       zint, 'r-.', 'Parent', ax(2)); 
         h7 = plot(lcwat_crm1D_m(iil,:), zint, 'r-.', 'linewidth',2, 'Parent',ax(2));
         end
         % plot dtcond spdt
         %
         %%%%%%%%%%%%%%%%%%%%%%%
         %   LS HEATING        %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_heating)
         ax(3) = axes('Position', ax1p,'YTickLabel','','XAxisLocation','bottom',...
                      'Xcolor','k','Color','none','Xlim',[-100 100],'Ylim',[zint(1) zint(end)],...
                      'XTick', xtickval3); hold on
%         h5 = plot(diffdt_m(iil,:,1), zint, 'Color', 'k', 'Parent', ax(3)); 
         h5 = plot(dtcond_cam1D_m(iil,:,1), zint, 'Color', 'k', 'Parent', ax(3)); 
         h6 = plot(spdt_cam1D_m(iil,:,1),   zint, 'Color', 'k', 'linewidth',1.5, 'Parent',ax(3)); 
         end
         %
         %%%%%%%%%%%%%%%%%%%%%%%
         %   STABILITY         %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_stability)
         ax(4) = axes('Position', ax1p, 'Color','none', 'XTickLabel', '','YTickLabel','', 'Xlim',[1 20],'Ylim',[0 16]); hold on
%         h7 = bar(1:2,cape(iil,1:2)/100,'Facecolor','none','EdgeColor','r','Parent',ax(4)); %xlim([1 20])
%         h8 = bar(3:4,cin(iil,1:2),'Facecolor','none','EdgeColor','b','Parent',ax(4));%xlim([1 20])
         h9 = bar(5:6,LI(iil,1:2),'Facecolor','none','EdgeColor','y','Parent',ax(4));%xlim([1 20])
         end
         %
         %%%%%%%%%%%%%%%%%%%%%%%
         %   RELHUM            %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_relhum)
         ax(5) = axes('Position',ax1p,'YTickLabel','','XAxisLocation','top',...
                      'Xcolor','b','Color','none','Xlim',[0 100],'Ylim',[zint(1) zint(end)],...
                      'XTick',xtickval4); hold on % <------- Make sure this hold on is here, otherwise it'll be overwritten
         h10= plot(rh_crm1D_m(iil,:), zint, 'b', 'linewidth', 2, 'Parent', ax(5));
         h11= plot(relhum_cam1D_m(iil,:), zint, 'b', 'linewidth', 1, 'Parent', ax(5));
         end
         %%%%%%%%%%%%%%%%%%%%%%%
         %   CRM RAIN          %
         %%%%%%%%%%%%%%%%%%%%%%%
         if (plot_crmrain)
         ax(6) = axes('Position', ax1p, 'Color','none', 'XTickLabel', '','YTickLabel','', 'Xlim',[1 32],'Ylim',[0 150]); hold on
         h13 = plot(rainrate(iil,:),'r-','Parent',ax(6));
         ax(7) = axes('Position', ax1p, 'Color','none', 'XTickLabel', '','YTickLabel','', 'Xlim',[1 32],'Ylim',[0 3]); hold on
         h12 = plot(cloudtype(iil,:),'*','Parent',ax(7));
         end
         %plot(spdt,zint,'Parent',ax2,'Color','g'); xlim([0 100]);
         %streamline(stream2(X,Z,squeeze(mean(var(iil,:,:,iit,iillt,1),5))',squeeze(mean(var(iil,:,:,iit,iillt,2),5))',X,Z)); % last two arg are startx startz
%         u_zave =
%         zprofile_centralize(squeeze(mean(mean(var(iil,:,:,iit,iillt,1),2),5)),parm.nx,parm.nzi);
%         w_zave = sparse(parm.nx,parm.nzi);
%         quiver(X,Z,squeeze(mean(var(iil,:,:,iit,iillt,1),5))',squeeze(mean(var(iil,:,:,iit,iillt,2),5))',0.3); hold off
%          quiver(X,Z,u_zave',w_zave',3); hold off; % this is for MATLAB2015
         % quiver(X,Z,u_zave',w_zave',3,'HeadStyle','plain'); hold off; pause % this is for MATLAB2016b
         iii=iii+1;%pause
      end; 
   end;
%   figure;plot(squeeze(mean(mean(stabilityIndex(:,:,iillt,idv.stabilityIndex.du3km_crm),1),3))*1000)
%   pause   
%   set(ha(1:end),'XTickLabel','')
%   set(ha(setdiff([1:nlons*nt],[1:nlons:nlons*nt])),'YTickLabel',''); % remove all Ytick except the left side subplots
%   print(gcf,['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/18km_lons11_lead0-10-lead48/organizeQcQiQr_basin' num2str(ib) '_t' num2str(it)],'-djpeg99')
%   savefig(gcf,['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/18km_lons11_lead0-10-lead48/organizeQcQiQr_basin' num2str(ib) '_t' num2str(it) '.fig'])
   toc
end % ib

function var = kmeans_kmedoids3D(iv, scale, var, ib, iit, iillt, nsamp, parm, kc, kmean_medoid)
   evalc(sprintf('[~, var] = %s(reshape(squeeze(var(:,:,:,iit,iillt,iv))*scale, parm.nlons*parm.nx*parm.nzi, nsamp)'',%d)', kmean_medoid, kc));
   var = reshape(permute(var,[2 1]),parm.nlons, parm.nx, parm.nzi, kc);

function var = kmeans_kmedoids2D(iv, scale, var, ib, iit, iillt, nsamp, parm, kc, kmean_medoid, maxiter)
   %evalc(sprintf('[~, var] = %s(reshape(squeeze(var(:,:,iit,:,iv))*scale, parm.nlons*parm.nzi, nsamp)'',%d)', kmean_medoid, kc));
   evalc(sprintf('[~, var] = %s(reshape(squeeze(var(:,:,iit,iillt,iv))*scale, parm.nlons*parm.nzi, nsamp)'',%d,''MaxIter'',maxiter)', kmean_medoid, kc));
   var = reshape(permute(var,[2 1]),parm.nlons, parm.nzi, kc);

function var = mean_median2D(iv, scale, var, ib, iit, iillt, mean_median) 
% Purpose: average/median over space-time index for a 2D CRM domain
% indices:: iv: variable, iillt: space-time, iit: lag-time, ib: basin 
   evalc(sprintf('var  = squeeze(%s(var(:,:,:,iit,iillt,iv),5))*scale',mean_median));

function var = mean_median1D(iv, scale, var, ib, iit, iillt, mean_median)
   if (isempty(iv))
      evalc(sprintf('var = squeeze(%s(var,2))',mean_median)); % zonal averaging
   elseif (ndims(var)==5)
      evalc(sprintf('var  = squeeze(%s(var(:,:,iit,iillt,iv),4))*scale',mean_median));
   end
function var = mean_median0D(iv, scale, var, ib, iit, iillt, mean_median)
% Purpose: average/median over space-time index iillt for a point
% variable, e.g., CAPE, CIN, LI
   evalc(sprintf('var  = squeeze(%s(var(:,iit,iillt,iv),3))*scale',mean_median)); 


function var = diffmean_median1D(iv, iv2, scale, var, ib, iit, iillt, mean_median)
   evalc(sprintf('var  = squeeze(%s(var(:,:,iit,iillt,iv)-var(:,:,iit,iillt,iv2),4))*scale',mean_median));

function  var = zprofile_centralize(varin,nx,nz) 
% centralize a vertical heating profile on a field map
   ix = ceil(nx/2); % the center point to 
   var = sparse(nx,nz);
   var(ix,:) = varin;


%function [ WRF_DBZ, WRF_DBZ_MAX ] = wrf_dbz(PRES,TK,QRAIN,QGRAUP,QSNOW,QVAPOR)
% calc radar reflectivity for precip rate estimation
% https://www.ncl.ucar.edu/Support/talk_archives/2012/att-0941/wrf_user_dbz.f   

function [T Z] = myclusterdata(X,distance,method,criterion,cutoff)
Y = pdist(X,distance); 
% distance : 
% 'euclidean', 'squaredeuclidean', 'seuclidean', 
% 'cityblock', 'minkowski', 'chebychev', 
% 'mahalanobis', 'cosine', 'correlation', 
% 'spearman', 'hamming', 'jaccard'
Z = linkage(Y,method); % k-row = k-link := horiz lines in dendrogram
% method : 
% 'average', 'centroid','complete','median', 
% 'single', 'ward', 'weighted'
T = cluster(Z,'criterion',criterion,'cutoff',cutoff);
% criterion : 
% 'inconsistent' or 'distance'

