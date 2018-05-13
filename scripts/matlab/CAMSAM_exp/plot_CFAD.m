function [svec variance] = plot_CFAD(vname,stats,ilndocn,ilatzone,season,dlat,plotmean,convstrat,signflip,dosave)
% vname='spdt';season='JJA';stats='mean';dlat='lat9090';convstrat=2;dosave=0;ilatzone=0;for ilndocn=[0 1]; plot_CFAD(vname,stats,ilndocn,ilatzone,season,dlat,0,convstrat,[],dosave);end
%function [pdf,cdf] = plot_CFAD(vname,season,stats,ilndocn,convstrat,dosave)
if ~exist('dosave')
  dosave=0;
end
%stats='mean'; 
thrhld = 25;
%dothr = input(['do variable conditioned on rainrate threshold = ' num2str(thrhld) '? (1/0)'])
dothr = 0;
%plotmean = input(['plot the ' stats '? (1/0)']);
%plotmean=1;
lndocnstr={'ocn','lnd','lndocn'};
lndocnStr={'Ocean','Land','Coast'};
startendstr = {'0','mid','end'};
pctlifestr = {'0%','50%','100%'};
ilndocn = ilndocn + 1;
%casei = 'F_2000_SPCAM_m2005_3hrly_f09f09_branch'
casei = 'F_2000_SPCAM_m2005_3hrly2';
latzone = {'trop','midlat','pole'};
%dlat = 'lat2525';
diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season]
diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/' lndocnstr{ilndocn} '/' latzone{ilatzone+1} '/']
%diri = ['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/f09f09/'];
%diro = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season '/' lndocnstr{ilndocn} '/'];
dircbar = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/'];
%diri = '/projects/rsmas/kirtman/gchen/cesm_spcam/archive/F_2000_SPCAM_m2005_3hrly1/atm/hist/'; 
load([diri '/mcs_cluster_var/mcs_cluster_' vname]);
load([diri '/mcs_clusters_1.mat'],'latC','nCl'); % lat zone: trop(0),midlat(1),polar(2)
eval(sprintf('v = %s;',vname));
load([diri '/mcs_cluster_stats/mcs_stats_' lndocnstr{ilndocn}],'iclndocn','nlndocnCl','iclndocn2');
iC = intersect(iclndocn2',find(latC==ilatzone)); % cluster indices conditioned on lnd/ocn and lat zone
load([diri '/mcs_cluster_parm.mat'])
load([diri '/mcs_cluster_var/mcs_cluster_ctype'])
load([diri '/mcs_cluster_var/mcs_cluster_prate'])
load([diri '/mcs_cluster_var/mcs_cluster_mcsllx4Cl'])
ipm = parm.nz; %maximum hybrid level
hp = flipud(parm.lev);
hp = hp(1:ipm);
% if ~ismember(convstrat,[0 1 2]) % specify prate instead of using ctype as criteria
%   pratethr = 50; % 16 > heavy rain > 4; 50 > very heavy rain > 16
binwidth=0.5;
   vrange={[-100 100],[-100 100],[-100 100],[-100 200],[-100 200],[-200 250],[-5 5],[30 100],[1e-5 1e-3],...
   [-100 100],[-100 100],[-100 100],[-100 100],[-100 100],...
   [-100 200],[-100 200],[-100 200],[-100 200],[-100 200],...
   };%0.01g/kg = 0.01e-3/kg = 
   xrange={[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-5 5],[30 100],[1e-5 1e-3],...
          [-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],[-100 100],};%0.01g/kg = 0.01e-3/kg = 
%{
   if convstrat==0
     vrange={[30 100],[30 100],[0 250],[0 250],[-750 250],[-50 100],[-5 5],[30 100],[1e-5 1e-3]};%0.01g/kg = 0.01e-3/kg = 
   elseif any(convstrat==[1,2])
     vrange={[-750 250],[-50 100],[-10 15],[30 100],[1e-5 1e-3]};%0.01g/kg = 0.01e-3/kg = 
   else
     vrange={[-750 250],[-50 100],[-8 8],[30 100],[1e-5 1e-3]};
   end
%}
%   xrange={[-150 150], [-50 100],[-8 8], [30 100], [1e-5 1e-3]};
%   vrange={[-8 15],[30 100],[1e-5 1e-3]};
   if dothr
     thr = ['g' num2str(thrhld) '_'];
   else
     thr = '';
   end 
   if convstrat == 0 & ~plotmean
     fileext = ['_precl_' thr];
   elseif convstrat == 1 & ~plotmean
     fileext = ['_precc_' thr];
   else % if convstrat == 2
     fileext = ['_precc_precl_' thr];
   end
% end
 tq = linspace(0,100,10);
 vnames = {'spdt','camdt','ddt_cam_sp','spdq','camdq','ddq_cam_sp','w_crm',...
           'rh_crm','qi_crm','zmdt','evaptzm','cmfdt','mpdt','macpdt','zmdq','evapqzm','cmfdq','mpdq','macpdq'};
% scale=[86400,86400,-86400,86400*3.6e3,86400*3.6e3,-86400*3.6e3,1,1,1];
 scale=[86400,86400,-86400,-86400*2260,-86400*2260,86400*2260,1,1,1,...
        86400,86400,86400,86400,86400,-86400*2260,-86400*2260,-86400*2260,...
        -86400*2260,-86400*2260];
 vdim = containers.Map(vnames,[1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1]);
 vmap = containers.Map(vnames,1:numel(vnames));
 iv = vmap(vname);
 cbarstr='%';
 caxisrange={[0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 10],[0 10],[0 10],...
             [0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 5],[0 5]};
 xlabels={'Q1c [K day^{-1}]','Q1c [K day^{-1}]','Q1c [K day^{-1}]','Q2 [K day^{-1}]','Q2 [K day^{-1}]','Q2 [K day^{-1}]','vertical velocity [m s^{-1}]','relative humidity [%]','ice mixing ratio [kg kg^{-1}]',...
         'Q1c [K day^{-1}]','Q1c [K day^{-1}]','Q1c [K day^{-1}]','Q1c [K day^{-1}]','Q1c [K day^{-1}]',...
         'Q2 [K day^{-1}]','Q2 [K day^{-1}]','Q2 [K day^{-1}]','Q2 [K day^{-1}]','Q2 [K day^{-1}]'};
 clev = linspace(caxisrange{iv}(1),caxisrange{iv}(2),100);
 nbin = diff(vrange{iv})/binwidth
% if convstrat == 0
%   nbin = 500;
% end
 v1{1}=[];v1{2}=[];v1{3}=[];
 for iic = 1:numel(iC)
   ic = iC(iic);
   if isempty(ctype{ic}); continue; end
   iit = [1,round(numel(ctype{ic})/2),numel(ctype{ic})]; % 0%,50%,100%
   for it = 1:numel(iit) 
     for ill=1:size(ctype{ic}{iit(it)},3)
       if vdim(vname)==2
         if ismember(convstrat,[0,1]) % condition on conv or strat
           if any(intersect(find(prate{ic}{iit(it)}(:,ill)>=thrhld),mcsllx4Cl{ic}{iit(it)}{ill})) & dothr
               v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(ctype{ic}{iit(it)}(:,ill)==convstrat),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
           elseif dothr==0
             v1{it} = cat(1,v1{it},v{ic}{1}(ctype{ic}{iit(it)}(:,ill)==convstrat,:,ill));
           end
         elseif convstrat ==2 % conv & strat rain mixed
           v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(~isnan(ctype{ic}{iit(it)}(:,ill))),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
         else % specify strong conv rainrate threshold
           v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(prate{ic}{iit(it)}(:,ill)>=convstrat),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
         end
       else
         v1{it} = cat(1,v1{it},v{ic}{iit(it)}(:,ill)'*scale(iv));
       end
     end 
   end
 end
  
%{   
 for nt2=3:numel(nlndocnCl); 
   if nlndocnCl{nt2}==0; continue; end
   for iic = 1:nlndocnCl{nt2};
     ic = iclndocn{nt2}(iic); 
     iit = [1,round(numel(ctype{ic})/2),numel(ctype{ic})]; % 0%,50%,100%
     for it = 1:numel(iit) 
       for ill=1:size(ctype{ic}{iit(it)},3)
         if vdim(vname)==2
           if ismember(convstrat,[0,1]) % condition on conv or strat
             if any(intersect(find(prate{ic}{iit(it)}(:,ill)>=thrhld),mcsllx4Cl{ic}{iit(it)}{ill})) & dothr
                 v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(ctype{ic}{iit(it)}(:,ill)==convstrat),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
             elseif dothr==0
               v1{it} = cat(1,v1{it},v{ic}{1}(ctype{ic}{iit(it)}(:,ill)==convstrat,:,ill));
             end
           elseif convstrat ==2 % conv & strat rain mixed
             v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(~isnan(ctype{ic}{iit(it)}(:,ill))),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
           else % specify strong conv rainrate threshold
             v1{it} = cat(1,v1{it},v{ic}{iit(it)}(intersect(find(prate{ic}{iit(it)}(:,ill)>=convstrat),mcsllx4Cl{ic}{iit(it)}{ill}),:,ill));
           end
         else
           v1{it} = cat(1,v1{it},v{ic}{iit(it)}(:,ill)'*scale(iv));
         end
       end 
     end
   end
 end
%   figure;
%}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the mean vertical profile for convective and stratiform case.
%{
if plotmean
%  for iz = 1:ipm
%    [h(iz) p(iz)] = ttest2(v1{1}(:,iz),v1{end}(:,iz),'Alpha',0.05);
%  end
  if convstrat ==0; dname = 'stratiform'; else; dname = 'convective'; end
  for it=1:numel(iit)
%    tmp = nan(2,ipm);
%    tmp(it,find(h)) = mean(v1{it}(:,find(h)),1);
    lw = 5; % linewidth
%    if convstrat == 0; lw = 3; hold on; else; lw = 3; figure('units','normalized','outerposition',[0 0 1 1]); end
    if convstrat == 0
      marker = {'b-.','b-','b:'};
    else
      switch vname
      case {'zmdt','zmdq'}
        marker = {'r-','r-.','r:'};
      case {'evaptzm','evapqzm'}
        marker = {'g-*','g-.','g:'};
      case {'mpdt','mpdq'}
        marker = {'b--','b-.','b:'};
      case {'macpdt','macpdq'}
        marker = {'y-.','y-.','y:'};
      case {'cmfdt','cmfdq'}
        marker = {'m-o','m-.','m:'};
      end
    end
    if strcmp(stats,'mean')
      [svec variance pc] = svdVar(v1{it}(:,1:ipm));
      variance = round(variance/sum(variance)*100);
%      [rvec varnew]   = rotate_pca(svec(:,1:3),pc(:,1:3));
%      varnew
      plot(svec(:,1)*signflip,hp,marker{it},'linewidth',lw,'DisplayName',dname); hold on;
%      plot(mean(v1{it}(:,1:ipm),1)./std(mean(v1{it}(:,1:ipm),1)),hp,marker{it},'linewidth',lw,'DisplayName',dname); hold on;
    elseif strcmp(stats,'median')
      plot(median(v1{it}(:,1:ipm),1)./std(mean(v1{it}(:,1:ipm),1)),hp,marker{it},'linewidth',lw,'DisplayName',dname); hold on;
    end
%    ha = plot(mean(v1{it}(:,find(h)),1),parm.zint(find(h))/1000,'r.','Markersize',lw*10);
%    ha = plot(zeros(1,ipm),hp,'k','linewidth',0.5); hold on
    ha = plot(zeros(1,ipm),hp,'k','linewidth',0.5); hold on
    set(get(get(ha,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%    plot(tmp(it,:),hp,'r','linewidth',lw,'r.','Markersize',15);
    xlabel(xlabels{iv})
%    ylabel('height [km]')
%    yticks([2:2:14])
%    xlim([-0.2 0.8])
%    ylim([0 15])
    ylabel('hybrid level','FontSize',16)
    set(gca,'YTick',linspace(200,1000,5));
    set(gca,'YTickLabel',linspace(200,1000,5));
    ylim([0 1000])
    set(gca,'Ydir','reverse');
    set(gca,'Fontsize',44)
return
    pause
  end
  title([season ' ' lndocnStr{ilndocn}])% ' (' pctlifestr{it} ' lifetime, n=' num2str(n(it)) ')']) 
%{
  dosave = 0;
    if dosave
      lg = legend('conv (0% life)','conv (50% life)','conv (100% life)','strat (0% life)','strat (50% life)','strat (100% life)')
      set(lg,'Fontsize',44)
      set(gcf,'color','w')
      box off
      legend(gca,'show'); 
      disp(['saveas(gcf,''' diro vname fileext stats '.fig'')'])
      disp(['saveas(gcf,' diro vname fileext stats '.png)'])
      disp(['crop(' diro vname fileext stats '.png)'])
    end
%}
% pause
end
%}




% vname = 'w_crm'; season='JJA'; ilndocn=1; convstrat=30; dosave=0; plot_CFAD(vname,season,ilndocn,convstrat,dosave)
% make sure to manually click figures to hold on the 0 and end figures
% vname = 'w_crm'; season='JJA'; ilndocn=1; convstrat=0; dosave=1; plot_CFAD(vname,season,ilndocn,convstrat,dosave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the contour CFAD 
 for it=1:numel(iit) % start and end of mcs lifecycle
   n(it) = size(v1{it},1);
   [pdf{it}, cdf{it}, vbincenters{it}] = binVarWithHeight(permute(v1{it}(:,1:ipm),[2,1]),vrange{iv},nbin,ipm);
%   pdf{it}(isnan(pdf{it})) = 0;
 end
% pvec = 0:5:100; % percentage vector
 for it=1:numel(iit)
   figure('units','normalized','outerposition',[0 0 1 1])
   hold on;
   contourf(vbincenters{it},hp,pdf{it}'*100,clev,'linestyle','none'); 
%   set(get(h,'title'),'string','[%]')
   colormap(cmap(12))
   caxis(caxisrange{iv})
  % for ip = 1:numel(pvec)-1 
  %   contour(vbincenters{it},hp,cdf{it}',[pvec(ip) pvec(ip)],'k');
  % end
   [~,h1] = contour(vbincenters{it},hp,cdf{it}'*100,[1 1],'k:','linewidth',1.5);
   contour(vbincenters{it},hp,cdf{it}'*100,[99 99],'k:','linewidth',1.5);
   [~,h2] = contour(vbincenters{it},hp,cdf{it}'*100,[10 10],'k:','linewidth',3.5);
   contour(vbincenters{it},hp,cdf{it}'*100,[90 90],'k:','linewidth',3.5);
 %  contour(vbincenters{i},hp,cdf{i}'*100,[25 25],'k','linewidth',1.5);
 %  contour(vbincenters{i},hp,cdf{i}'*100,[75 75],'k','linewidth',1.5);
 %  xticks(vtick{iv});
 %  yticks(hp)
    lw = 5; % linewidth
%    if convstrat == 0; lw = 3; hold on; else; lw = 3; figure('units','normalized','outerposition',[0 0 1 1]); end
  if convstrat ==0; dname = 'stratiform'; else; dname = 'convective'; end
  if convstrat == 0
    marker = {'b-.','b-','b:'};
  else
%    marker = {'r-.','r-','r:'};
    marker = {'k-','k-','k-'};
  end
  if strcmp(stats,'mean')
   h3= plot(mean(v1{it}(:,1:ipm),1),hp,marker{it},'linewidth',lw,'DisplayName',dname); hold on;
%  elseif strcmp(stats,'median')
   h4= plot(median(v1{it}(:,1:ipm),1),hp,marker{it},'linewidth',lw,'color',[.5 .5 .5],'DisplayName',dname); hold on;
  end
  if ilndocn==2 & ismember(vname,{'spdt','spdq'})
    hl = legend([h4,h3,h2,h1],'median',stats,'10 & 90 prctile','1 & 99 prctile');
    set(hl,'Box','off') 
  end
  ha = plot(zeros(1,ipm),hp,'k','linewidth',0.5); hold on
  title([season ' ' lndocnStr{ilndocn} ' (' pctlifestr{it} ' lifetime)']) 
%  ylabel('height [km]','FontSize',16)
  ylabel('hybrid level','FontSize',16)
  set(gca,'YTick',linspace(200,1000,5));
  set(gca,'YTickLabel',linspace(200,1000,5));
  ylim([0 1000])
  set(gca,'Ydir','reverse')
  xlabel(xlabels{iv},'FontSize',16)
%  yticks([2:2:14])
  if convstrat==0 
    xlim([-3 3])
  else
    xlim([-5 8])
  end
  xlim(xrange{iv})
  set(gca,'Layer','top')
  set(gca,'fontsize',44)  
  set(gcf,'color','w')
  if dosave
%    saveas(gcf, [diro 'CFAD_' vname fileext  startendstr{it} '.fig'])
    saveas(gcf, [diro 'CFAD_' vname fileext  startendstr{it} '.png'])
    crop([diro 'CFAD_' vname fileext  startendstr{it} '.png'])
  else
    disp(['saveas(gcf,''' diro 'CFAD_' vname fileext  startendstr{it} '.png'')'])
    disp(['crop(''' diro 'CFAD_' vname fileext  startendstr{it} '.png'')'])
  end
  docbar = 0;
  if docbar
    figure('units','normalized','outerposition',[0 0 1 1])
    h = colorbar('southoutside'); 
    axis off
    colormap(cmap(12))
    caxis(caxisrange{iv})
    set(gca,'fontsize',44);
    h.Label.String = '[%]';
%    h.Label.Position = [195 0 0];
  h.Label.Position = [caxisrange{iv}(2)+0.08*diff(caxisrange{iv}) 0 0];
    dosave = 1;
    if dosave 
      saveas(gca,[dircbar '/CFAD_colorbar.fig'])
      saveas(gca,[dircbar '/CFAD_colorbar.png'])
      crop([dircbar '/CFAD_colorbar.png'])
    end
  end
 end
 %eval(sprintf('fo = ''CFAD_%s.mat'';',vname))
 %save(fo,'pdf','cdf','vbincenters');

% vname='w_crm';season='JJA';stats='mean'; ilndocn=0;convstrat=0;dosave=0;plot_CFAD(vname,season,stats,ilndocn,convstrat,dosave)
%  vname='spdt';season='JJA';stats='median';convstrat=2;dosave=1; for ilndocn=[0 1]; plot_CFAD(vname,season,stats,ilndocn, convstrat,dosave,0,0);end


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Purpose: Mean plot of CAM ZMDT,EVAPTZM,CMFDT,MPDT,MACPDT
%tq='t';
tq='q';
figure('units','normalized','outerposition',[0 0 1 1]);
vname = {['zmd' tq],['evap' tq 'zm'],['cmfd' tq],['mpd' tq],['macpd' tq]};
%{
if strcmp(tq,'q')
  signflip = [1,-1,1,1,1];
elseif strcmp(tq,'t')
  signflip = [-1,-1,1,1,1];
end
%}
season='JJA';stats='median'; ilndocn=0;convstrat=2;
casei = 'F_2000_SPCAM_m2005_3hrly2';
dlat = 'lat2525';
lndocnstr={'ocn','lnd','lndocn'};
diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/' lndocnstr{ilndocn+1} '/']
%for iv=1:numel(vname); [svec{iv},var{iv}] =plot_CFAD(vname{iv},season,stats,ilndocn,convstrat,0,1,signflip(iv)); end
for iv=1:numel(vname); plot_CFAD(vname{iv},season,stats,ilndocn,convstrat,0,1); end
if strcmp(tq,'q')
  xlabel('normalized Q2')
%  xlabel('EOF1 Q2')
elseif strcmp(tq,'t')
  xlabel('normalized Q1c')
%  xlabel('EOF1 Q1c')
end
xlim([-4 4])
%xlim([-1 1])
ha = legend('P_{deep}','P_{evap}','P_{shallow}','P_{micro}','P_{macro}')
%{
ha = legend(['P_{deep}    (%Var=' num2str(var{1}(1)) ')'],...
            ['P_{evap}    (%Var=' num2str(var{2}(1)) ')'],... 
            ['P_{shallow} (%Var=' num2str(var{3}(1)) ')'],... 
            ['P_{micro}   (%Var=' num2str(var{4}(1)) ')'],... 
            ['P_{macro}  (%Var=' num2str(var{5}(1)) ')'])
%}
set(ha,'box','off')
set(ha,'location','northwest')
saveas(gcf, [diro 'camd' tq '_decomp.png'])
crop([diro 'camd' tq '_decomp.png'])
%saveas(gcf, [diro 'camd' tq '_eofdecomp.png'])
%crop([diro 'camd' tq '_eofdecomp.png'])





%}
