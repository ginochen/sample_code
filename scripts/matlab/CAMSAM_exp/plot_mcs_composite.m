function plot_mcs_composite(vname,stats,ilndocn,ilatzone,season,dlat,imode,dosave)
%function plot_mcs_composite(vname,stats,ilndocn,season,imode,dosave)
%  figure('units','normalized','outerposition',[0 0 1 1])
% ddt ddq calc according to ~/notes/matlab/cellfun.txt
  if ~exist('dosave')
    dosave=0
  end
  doqtile=0; % plot 25 and 75 quartile
  docbar = 0;
  doyy = 0;  
%  figure('units','normalized','outerposition',[0 0 1 1])
  casei='F_2000_SPCAM_m2005_3hrly2'
%  casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch'
%  dlat = 'lat2525';
  diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season]
  res = 'f09f09';
  latzone = {'trop','midlat','pole'};
  vnames0 = {'prect','precc','precl','freq_precc','freq_precl',... % prect = precc_cam
             'dU03km_cam','du03km_abs','dv03km_abs','frac_precc','frac_precl','frac_prec0','frac_prec3','frac_prec5','enst','flnt','B','Bave','LIave'}; %16
  vnames1 = {'spdt','spdq','spmc','spmcup','spmcdn','div',... % 6
             'w_strat','w_conv','qv_strat','qv_conv','qi_strat','qi_conv',... % 12
             'qc_strat','qc_conv','qr_strat','qr_conv','rh_strat','rh_conv',... %18 
             'u_strat','u_conv','camdt','camdq','ddt_cam_sp','ddq_cam_sp'}; %22
  vdim = containers.Map([vnames0,vnames1],[zeros(1,numel(vnames0)),ones(1,numel(vnames1))]);
  ivv  = containers.Map([vnames0,vnames1],[1:numel(vnames0),1:numel(vnames1)]);
  statname = {'mean','median','std','pca','rpca'};
  istats = containers.Map(statname,[1:numel(statname)]);
  variancethr = 0.6; % for pca truncated
  maxnt = 20;% 7;
%  ikm = 20; %17; % 10km height index is 12
%  ikm = 17; %17; % 10km height index is 12
  L = 2260; % latent heat of vaporization
  lndocnstr = {'ocn','lnd','lndocn','nocond'};
  lndocnStr = {'Ocean','Land','Coast'};
  ilndocn = ilndocn + 1; % +1 for actual index
  load([diri '/mcs_cluster_parm.mat'])
  load([diri '/mcs_clusters_1.mat'],'latC','nCl'); % lat zone: trop(0),midlat(1),polar(2)
  ikm = parm.nz; 
  hp = flipud(parm.lev);
  hp = hp(1:ikm); 
  diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/' lndocnstr{ilndocn} '/' latzone{ilatzone+1} '/'] % output dir for figures
  if vdim(vname) == 0
%    vname = {'precc','precl','freq_precc','freq_precl','du03km_abs','dv03km_abs'};
    ylimrange = {[0 5],[0 5],[0 5],[0 .8],[0 .8],[0 15],[0 9.5],[0 10],[0 1],[0 1],[0 1],[0 1],[0 1],[0 0.003],[-0.1 0.1],[40 100],[40 100],[-7 5]};
    yticks    = {[],[],[],linspace(0,.8,5),linspace(0,.8,5),round([linspace(0,9.5,5)],1),round([linspace(0,9.5,5)],1),linspace(0,12.5,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6),linspace(0,0.003,6),'',linspace(-0.1,0.1,6),40:10:100,40:10:100,-7:1:1};
%    ylimrange = {[0 5],[0 5],[0 .8],[0 .8],[0 45],[0 10],[0 1],[0 1],[0 1],[0 1],[0 1]};
%    yticks    = {[],[],linspace(0,.8,5),linspace(0,.8,5),[linspace(0,12.5,6),45],linspace(0,12.5,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6),linspace(0,1,6)};
    scale = [1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
    if ilndocn == 1 %ocn
      mark  = {'-','-','-','o','o','-','-','--','-','-','-','-','-','-','-','-','-','-','-'}; %16
    else
      mark  = {'-','-','-','o','o','-','-','--','-.','-.','-.','-.','-.','-','-','-','-','-','-'};
    end
    color = {'c','b','c','k','k','r','r','b','r','b','b','b','r','k','k','k','k','k','k'};
    ylabels = {'','','','rain fraction','rain fraction','|\DeltaU| [m s^{-1}]','|\DeltaU| [m s^{-1}]','|\DeltaU| [m s^{-1}]','fractional rainrate [%]','fractional rainrate [%]','fractional rainrate [%]','fractional rainrate [%]','fractional rainrate [%]','m^2 s^{-2}','','m s^{-1}','m s^{-1}','K'};
  elseif vdim(vname) == 1
%    vname = {'spdt','spdq','spmc','spmcup','spmcdn','div'}; 
%    ylimrange = [0 parm.lev(ikm)]
  
    scale = [86400, -2260*86400, 1, 1, 1, 1,1,1,1,1e3,1e3,1e3,1e3,1e3,1e3,1e3,1,1,1,1,86400,-2260*86400,-86400,2260*86400]; %22
%    scale = [86400, -L*86400, 1, 1, 1, 1,1];
    cbarunits = {'[K day^{-1}]','[K day^{-1}]','[kg m^{-2} s^{-1}]','[kg m^{-2} s^{-1}]','[kg m^{-2} s^{-1}]','[s^{-1}]',... %6
                 '[m s^{-1}]','[m s^{-1}]','[g kg^{-1}]', '[g kg^{-1}]','[g kg^{-1}]', '[g kg^{-1}]',... %12
                 '[g kg^{-1}]', '[g kg^{-1}]','[g kg^{-1}]', '[g kg^{-1}]','[%]', '[%]',... %18
                 '[m s^{-1}]','[m s^{-1}]','[K day^{-1}]','[K day^{-1}]','[K day^{-1}]','[K day^{-1}]'}; %22
    cint = {[-20:2.5:40],[-20:2.5:40],linspace(-0.01,0.03,300),linspace(-0.01,0.03,300),linspace(0,0.03,300),linspace(-9e-7,9e-7,300),... %6
            [-0.25:0.025:0.25],[-0.25:0.025:0.25],linspace(0,0.03,300),linspace(0,0.03,300),linspace(0,7e-2,10),linspace(0,7e-2,10),... %12
            linspace(0,8e-2,10),linspace(0,1.3e-1,10),linspace(0,2.5e-1,10),linspace(0,2.5e-1,10),linspace(0,90,10),linspace(0,90,10)... %18
            linspace(-15,15,10),linspace(-5,5,10),[-20:2.5:40],[-20:2.5:40],[-20:2.5:40],[-20:2.5:40]}; %22
%    ncmap = [23,16,1,1,1,1,1,1,12,12,12,12,12,12,12,12,12]; 
    ncmap = [23,23,1,1,1,1,1,1,1,12,19,19,19,19,19,19,19,19,1,1,23,23,23,23]; %22
    flipm = [0,  0,0,0,0,0,0,0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0, 0, 0, 0,0];
    caxislim = {[-10 35],[-10 35],[-0.01 0.03],[-0.01 0.03],[0 0.03],[-8e-8 8e-8],... %6
                [-0.25 0.25],[-0.25 0.25],[0 3e-2],[0 3e-2],[0 7e-2],[0 7e-2],... %12
                [0 6e-2],[0 1.3e-1],[0 2.5e-1],[0 2.5e-1],[0 1e-1],[0 5e-1]... %18 
                [-15 15],[-5 5],[-10 35],[-10 35],[-10 35],[-10 35]}; %22
    color = containers.Map({'qc_strat','qi_conv','qc_conv'},{[0 0 .6],[.6 0 0],[0 0 .6]});
  end
  eval(sprintf('load(''%s/mcs_cluster_stats/mcs_stats_%s.mat'',''%s'',''nlndocnCl'',''ntvec2'',''iclndocn'',''iclndocn2'');',diri,lndocnstr{ilndocn},vname));
  if isempty(ilatzone)
    ilatC = 1:numel(iclndocn2); % entire cluster over land without lat conditions
  else
    ilatC = find(latC(iclndocn2)==ilatzone); %clusters in the tropics(0) midlat(1)
    % latC: lat zone condition for all clusters, iclndocn2: cluster index
    % that are over lnd/ocn, ilatzone: lat zone condition
  end
  eval(sprintf('vtmp = %s.mcscomposite;',vname));
  timerange = linspace(0,100,maxnt);
  if strcmp(stats,'pca')
    varitmp = vtmp.pcvariance;
    varitmp1 = sum(varitmp(1)/sum(varitmp));
    nm=1;
    while varitmp1 < variancethr % set truncation modes according to total variance ratio 
      nm = nm + 1;
      varitmp1 = sum(varitmp(1:nm)/sum(varitmp));
    end
    disp(['maximum mode when variance reaches ' num2str(variancethr) ' is ' num2str(nm) ]) 
    if ~exist('imode')
      disp('do truncated pca summation til the fractional variance reaches')
      pc_meantmp = vtmp.pc_mean;
      svectmp = vtmp.svec;
      truncm = nm; 
      sumvar = 0; 
      for im = 1:truncm 
        sumvar = sumvar + pc_meantmp(im)*svectmp(:,im); % summing the truncated variable
      end 
      vtmp=reshape(sumvar,parm.nzi,maxnt)*scale(ivv(vname)); 
    else
      vtmp = reshape(vtmp.svec(:,imode),ikm,maxnt);
      if vdim(vname) == 1
%        for i = 1:numel(imode)
%          vtmp1 = vtmp1 + reshape(scale(ivv(vname))*vtmp.pc_mean(imode(i))*vtmp.svec(:,imode(i)),ikm,maxnt);
%        end
        if ismember(vname,{'w_strat','w_conv','qc_strat','qc_conv','qi_strat','qi_conv','qr_strat','qr_conv'})
          contour(timerange,hp,vtmp,'k',round(cint{ivv(vname)}),2);
        else
%          ylim(ylimrange) 
          contourf(timerange,hp,vtmp,100,'linestyle','none');
          h= colorbar;
          set(get(h,'title'),'string',cbarunits(ivv(vname))); % units on top
          colormap(cmap(ncmap(ivv(vname)),flipm(ivv(vname)),18));
%          caxis([-0.17 0.7])
%          caxis([-0.01 0.01])
%          caxis(caxislim{ivv(vname)});
%          colormap(cmap(1))
          caxis([-.1 .1]) 
        end
%{
        title([season ' ' lndocnStr{ilndocn} ' (%Var = ' num2str(round(varitmp(imode)'/sum(varitmp)*100)) ', n = ' num2str(sum([nlndocnCl{:}])) ')'])
%}
      elseif vdim(vname) == 0
        mark = {'k-','k--','k:'};
        [~,loads] = sign_flip({vtmp.pc(:,1:imode),vtmp.svec(:,1:imode)},reshape(vtmp.data,maxnt,numel(vtmp.data)/(maxnt))');
        plot(timerange,loads{2}(:,imode),mark{ilndocn},'linewidth',3);
%        plot(timerange,signs(ilndocn)*vtmp.pc_mean(imode)*vtmp.svec(:,imode),mark{ilndocn},'linewidth',3);
      end      
      disp([num2str(imode) ' mode has variance portion of ' num2str(varitmp(imode)'/sum(varitmp))]); 
    end
%    ylabel('height [km]')
    ylabel('hybrid level')
%{
    xlabel('lifetime [%]')
%}
  elseif strcmp(stats,'rpca')
    [~,loads] = sign_flip({vtmp.pc(:,1:imode),vtmp.svec(:,1:imode)},reshape(vtmp.data,ikm*maxnt,numel(vtmp.data)/(ikm*maxnt))');
    [rvec varnew] = rotate_pca(loads{2}(:,1:imode),loads{1}(:,1:imode)); % (svec, pc)
    disp(['fractional variance for rotated pca modes = ' num2str(varnew/sum(varnew)*100)])
%    disp(['fractional variance for unrotated pca modes = ' num2str(vtmp.pcvariance(1:imode)'/sum(vtmp.pcvariance)*100)])
    for im=1:imode
      figure(im)
      set(im,'units','normalized','outerposition',[0 0 1 1]); 
%  subplot(1,2,1)
      contourf(timerange,hp,reshape(rvec(:,im),ikm,maxnt),100,'linestyle','none');
%  subplot(1,2,2)
%      contourf(timerange,hp,reshape(loads{2}(:,im),ikm,maxnt),100,'linestyle','none');
%      h = colorbar;
%      set(get(h,'title'),'string')); % units on top
      colormap(cmap(1));
      caxis([-0.1 0.1])
      set(gca,'FontSize',44)
      title([season ' ' lndocnStr{ilndocn} ' (REOF' num2str(im) ', %Var=' num2str(round((varnew(im)'/sum(varnew))*100,2)) ', n=' num2str(sum([nlndocnCl{:}])) ')'])
%      ylabel('height [km]')
      ylabel('hybrid level')
      xlabel('lifetime [%]')
%      figure;
%      contourf(timerange,hp,reshape(vtmp.svec(:,im),ikm,maxnt),100,'linestyle','none');
%      h = colorbar;
%      colormap(cmap(1));
%      caxis([-0.1 0.1])
    end
  else
    if vdim(vname)==1
      ind = {':',':',ilatC};
    elseif vdim(vname)==0
      ind = {':',ilatC};
    end 
      %eval(sprintf('vtmp = vtmp.%s*scale(ivv(vname));',stats));
      if strcmp(stats,'mean')
        vtmp = nanmean(vtmp.data(ind{:})*scale(ivv(vname)),vdim(vname)+2);
      elseif strcmp(stats,'median')
        vtmp = nanmedian(vtmp.data(ind{:})*scale(ivv(vname)),vdim(vname)+2);
      end 
      %vtmp = vtmp.data*scale(ivv(vname));
  end
  if any(strcmp(stats,{'mean','median'})) % stats excluding pca and rpca
    switch vdim(vname)
    case 0
      if doyy
        yyaxis right
        set(gca,'YTick',yticks{ivv(vname)},'YColor',color{ivv(vname)})
      else
        set(gca,'YTick',yticks{ivv(vname)})
      end
      if ilndocn==1 | ~strcmp(vname,{'du03km_abs'}); doslope=0; else; doslope=1; end % cal slope for low level shear evolution
      if strcmp(vname,'dU03km_cam'); vtmp = abs(vtmp); end
      plotquartile(timerange,vtmp',stats,mark{ivv(vname)},color{ivv(vname)},doqtile,doslope);
      xlabel('lifetime [%]')
      ylim(ylimrange{ivv(vname)})
      ylabel(ylabels{ivv(vname)})
      %title([ season ' ' lndocnStr{ilndocn} ' (n = ' num2str(sum([nlndocnCl{:}])) ')'])
%      title([ season ' ' lndocnStr{ilndocn} ' (n = ' num2str(numel(ilatC)) ')'])
    case 1
      if ismember(vname,{'w_strat','w_conv','qc_strat','qc_conv','qi_strat','qi_conv','qr_strat','qr_conv'})
        v = round(cint{ivv(vname)},2);
%        v1 = [v(1),v(3),v(5),v(7),v(8),v(9),v(10)]
        v1 = [v(1:2:end-5),v(end-4:end)]
        if ismember(vname,{'qi_conv','qi_strat','qc_conv','qc_strat'}) 
          clev=[0:0.02:0.08];
          for icl = 1:numel(clev) 
            contour(timerange,hp,vtmp(1:ikm,:),[clev(icl) clev(icl)],'color',color(vname),'linewidth',icl^1.5); hold on
          end
          return
        else
        [C, h] = contourf(timerange,hp,vtmp(1:ikm,:),cint{ivv(vname)},'linestyle','none'); hold on
        cin = cint{ivv(vname)};
        contour(timerange,hp,vtmp(1:ikm,:),[cin(cin<=0)],'b--');
        contour(timerange,hp,vtmp(1:ikm,:),[cin(cin>0)],'r-');
        %clabel(C,h,v1,'LabelSpacing',1500,'Fontsize',30)%,'BackgroundColor','w')
%        if ismember(vname,{'w_strat','w_conv'})
%          clabel(C,[v(1:10:numel(v)),v(numel(v)-1:numel(v))],'Fontsize',30)
        end
%        else
%          clabel(C,[v(1:5:numel(v)),v(numel(v)-1:numel(v))],'Fontsize',30,'Color','Red')
%        end
      else
        contourf(timerange,hp,vtmp(1:ikm,:),cint{ivv(vname)},'linestyle','none'); 
%      h= colorbar;
%      set(get(h,'title'),'string',cbarunits(ivv(vname))); % units on top
%      ylabel(h,cbarunits(ivv)); % units on the side
      end
%        vtmpn = vtmp(1:ikm,:);
%        vtmpn(vtmpn>0)=0; 
      hold on
      [~,hn]=  contour(timerange,hp,vtmp,[-20:2.5:0],'--','color',[0 0 0.6],'linewidth',0.5);hold on
%      hn.LevelStep = 1; % set contour interval
      colormap(cmap(ncmap(ivv(vname)),flipm(ivv(vname)),18));
      caxis(caxislim{ivv(vname)}); 
      set(gca,'YTick',linspace(200,1000,5));
      set(gca,'YTickLabel',linspace(200,1000,5));
      ylim([0 1000])
      set(gca,'Ydir','reverse')
%      set(gca,'YTickLabel',hp);
%      ylim(ylimrange) 
      %ylabel('height [km]')
      ylabel('hybrid level')
      xlabel('lifetime [%]')
    end
    %title([ season ' ' lndocnStr{ilndocn} ' (n = ' num2str(sum([nlndocnCl{:}])) ')'])
   % title([ season ' ' lndocnStr{ilndocn} ' (n = ' num2str(numel(ilatC)) ')'])
  end
  if strcmp(stats,'pca') 
%    if ilndocn==1
%      xlim([-1 1])
%{
      if ismember(vname,{'spdt','camdt','ddt_cam_sp'})
        xlim([-30 40])
      else
        xlim([-15 10])
      end
      set(gca,'XTick',[-30 -20 -10 -5 0 5 10 20 30 40])
%}
%    elseif ilndocn==2
%      xlim([-30 40])
%      xlim([-30 10])
%      set(gca,'XTick',[-30 -20 -10 -5 0 5 10 20 30 40])
%    end
  else
    set(gca,'XTick',[0:20:100])
  end
  set(gcf,'color','w')
  set(gca,'FontSize',44)
  set(gca,'Layer','top')
  if any(strcmp(vname,{'spdt'}))
    if ismember(stats,{'rpca'})
      for im = 1:imode
        fname{im} = ['rEOF' num2str(im) '_' vname];
      end
    elseif ismember(stats,{'pca'})
      for im = 1:imode
        fname{im} = ['EOF' num2str(im) '_' vname];
      end
    elseif ismember(stats,{'mean','median'})
      fname = ['complife_' vname '_' stats];
    end
  elseif ismember(vname,{'frac_precc','frac_precl'})
    fname = ['complife_frac_precc_precl_' stats];
    if ismember(vname,{'frac_precl'})
      lg=legend('conv (ocean)','strat  (ocean)','conv (land)','strat  (land)'); set(lg,'box','off','location','southeast');
    end
  elseif all(strcmp(vname,{'spdq','du03km_abs'}))
    fname = ['complife_spdq_du03km_' stats]
  elseif ismember(vname,{'frac_prec0','frac_prec3','frac_prec5'})
    if ismember(vname,{'frac_prec5'})
%      lg=legend('>0.1 mm hr^{-1}','>3 mm hr^{-1}','>5 mm hr^{-1}'); set(lg,'box','off','location','southeast')
      lg=legend('>3 mm hr^{-1} (ocean)','>5 mm hr^{-1} (ocean)','>3 mm hr^{-1} (land)','>5 mm hr^{-1} (land)'); set(lg,'box','off','location','southeast')
    end
    fname = ['complife_frac_prec5_' stats];
  elseif ismember(vname,{'spmc'})
    if ismember(stats,{'rpca'})
      fname = 'rEOF_spmc';
    end 
  %
  %elseif ismember(vname,'')
  else
    fname = ['complife_' vname '_' stats];
  end
  if dosave
    if ismember(stats,{'rpca'})
      for im=1:imode
        saveas(im,[diro '/' fname{im} '.png'])
        crop([diro '/' fname{im} '.png'])
      end
    else   
      saveas(gcf,[diro '/' fname '.png'])
      crop([diro '/' fname '.png'])
    end
%    close all
  end
  disp(['saveas(gcf,' ''''  diro '/' fname '.png' '''' ')'])
  disp(['crop(' '''' diro '/' fname '.png' '''' ')'])
% create colorbar only plots
%
if docbar & vdim(vname) 
  figure('units','normalized','outerposition',[0 0 1 1]); 
  axis off; 
  h=colorbar('southoutside'); 
  colormap(cmap(ncmap(ivv(vname)),flipm(ivv(vname)),18)); 
  caxis(caxislim{ivv(vname)}); 
  h.Position = [0.05 0.2170 0.7750 0.0181];
  set(gca,'Fontsize',44); 
  h.Label.String = cbarunits{ivv(vname)}
  h.Label.HorizontalAlignment = 'right';
  h.Label.Position = [caxislim{ivv(vname)}(2)+0.225*diff(caxislim{ivv(vname)}) 0 0];
  saveas(gcf,['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/complife_' vname '_colorbar.png']); 
        crop(['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/complife_' vname '_colorbar.png'])
end


%%%%%%%%%%%%%%%%%%%%%% FUCNTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qtile = plotquartile(timerange,data,stats,mark,color,doqtile,doslope)
  % data(clustsamp,var)
  if strcmp(stats,'median')
    y =  nanmedian(data,1);
    plot(timerange,y,mark,'color',color,'linewidth',4); hold on
  elseif strcmp(stats,'mean')
    y =  nanmean(data,1);
    plot(timerange,nanmean(data,1),mark,'color',color,'linewidth',4); hold on
  end
  if doslope % regression slope
    [b,~,~,~,stats] = regress(data(:),[ones(numel(data),1),reshape(repmat(timerange,size(data,1),1),numel(data),1)]);
    text(timerange(10),y(10)+0.4,['slope = ' num2str(round(b(2),3))],'Fontsize',44,'Color','Red')
  end
  % stats:[R^2-stat,F-stat,F-stat p-val,error varaiance]
  % this is same as performing fitlm()
  % fitlm([reshape(repmat(timerange,size(data,1),1),numel(data),1)],data(:))

  if doqtile
    [qtile] = prctile(data,[25 75]);
    plot(timerange,qtile,'--','color',color,'linewidth',2)
  end

%{
dlat='lat9090'; season='JJA'; stats='median'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); hold on; plot_mcs_composite('B',stats,ilnd,season,dlat); end
season=''; stats='median'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); hold on; plot_mcs_composite('div',stats,ilnd,season); end


vname='spdt';season='JJA';ilatzone=0;dlat='lat9090';stats='mean'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season,dlat,[],dosave); end

vname='w_strat';season='';dlat='lat9090';stats='median'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season,imode,dosave); end

vname='w_conv';season='';dlat='lat9090';stats='median'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season,imode,dosave); end

vname='du03km_abs';season='';dlat='lat9090';stats='median'; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season,0,dosave); end

season=''; stats='median'; figure('units','normalized','outerposition',[0 0 1 1]);for ilnd=0:1; hold on; plot_mcs_composite('frac_prec0',stats,ilnd,season,imode);    plot_mcs_composite('frac_prec3',stats,ilnd,season); plot_mcs_composite('frac_prec5',stats,ilnd,season,0); end; title(season)

vname={'frac_precc','frac_precl'};season='JJA'; ilatzone=0; stats='mean'; dosave=0; figure('units','normalized','outerposition',[0 0 1 1]); for ilnd=0:1; hold on;  plot_mcs_composite(vname{1},stats,ilnd,ilatzone,season,dlat,[],dosave);  plot_mcs_composite(vname{2},stats,ilnd,ilatzone,season,dlat,[],dosave); end; 

vname='spdt'; stats='rpca'; season=''; imode=3; dosave=0; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season,imode,dosave);end

season=''; stats='median'; dosave=1; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); hold on; plot_mcs_composite('spdq',stats,ilnd,ilatzone,season); plot_mcs_composite('du03km_abs',stats,ilnd,ilatzone,season,[],dosave); end

vname='spdq';season=''; stats='median'; dosave=1; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); plot_mcs_composite(vname,stats,ilnd,ilatzone,season); end

season=''; stats=''; dosave=0; figure('units','normalized','outerposition',[0 0 1 1]); for ilnd=0:1; hold on; plot_mcs_composite('frac_prec0',stats,ilnd,ilatzone,season,imode); plot_mcs_composite('frac_prec3',stats,ilnd,ilatzone,season); plot_mcs_composite('frac_prec5',stats,ilnd,ilatzone,season,0,dosave); end; 

 imode=[2]; ilnd=0; season='JJA'; plot_mcs_composite('spdt','rpca',ilnd,ilatzone,season,imode);
 imode=[2]; ilnd=0; season='JJA'; plot_mcs_composite('spmc','rpca',ilnd,ilatzone,season,imode);

% load JJA/mcs_cluster_stats/mcs_stats_lnd spmc
% im = 3; contourf(reshape(spmc.mcscomposite.pc_mean(im)*spmc.mcscomposite.svec(:,im),17,20))

% Water phase
season={'DJF','JJA'};for is=1:2; stats='median'; vname='qc_strat'; dosave=1; for ilnd=0:1; figure('units','normalized','outerposition',[0 0 1 1]); hold on; plot_mcs_composite(vname,stats,ilnd,ilatzone,season{is},imode,  dosave); end; end

% for w_conv + qi_conv and w_strat + qc_strat
vname={'w_strat','qc_strat'};ilatzone=1;ilnd=0;season='JJA';figure('units','normalized','outerposition',[0 0 1 1]); stats='mean'; dosave=0; for iv=1:2;plot_mcs_composite(vname{iv},stats,ilnd,ilatzone,season,dlat,[]); hold on;end

vname={'w_conv','qc_conv'};ilatzone=1;ilnd=0;season='JJA';figure('units','normalized','outerposition',[0 0 1 1]); stats='mean'; dosave=0; for iv=1:numel(vname); plot_mcs_composite(vname{iv},stats,ilnd,ilatzone,season,dlat,[]); hold on;end

vname='qc_conv'; plot_mcs_composite(vname,stats,ilnd,ilatzone,season,0,dosave);

% plot shear EOF1
% season = 'DJF';
% season = 'JJA';
 figure('units','normalized','outerposition',[0 0 1 1]); 
 ilnd=0;hold on; imode = [1]; plot_mcs_composite('du03km_abs','pca',ilnd,season,ilatzone,imode); % make sure to press ctrl+c after each curve
 ilnd=1;hold on; imode = [1]; plot_mcs_composite('du03km_abs','pca',ilnd,season,ilatzone,imode);
% legend('ocean (%Var = 93)','land    (%Var = 84)') % for JJA
% legend('ocean (%Var = 95)','land    (%Var = 87)') % for DJF
 set(gca,'XTick',[0:20:100])
 set(gca,'Fontsize',44)
 xlabel('lifetime [%]')
 title([season]) 
 ylabel('|\DeltaU| EOF1')
% saveas(gcf,['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season '/EOF1_du03km.fig'])
 saveas(gcf,['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season '/EOF1_du03km.png'])
 crop(['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season '/EOF1_du03km.png'])
%}
