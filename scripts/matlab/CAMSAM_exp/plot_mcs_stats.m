function [ha vtmp1] = plot_mcs_stats(ivv,istat,ilndocn,ha,vtmp1)
%function [ha vtmp1] = plot_mcs_stats(ivv,istat,ilndocn,ha,vtmp1)
% ha : subplot handles, output to plot multiple variables on the same subplots
% vtmp1: output precc to use bar([vtmp,vtmp2],'stacked')
%%%%%%  PARAMETER  %%%%%%%%%%%%%%%%%%%
%  ivv: 1 spdt 2 dtcond 3 spdq 4 prect 5 precc 6 precl 7 divmax_crm 8 LI_crm
%  istat: 1 mean, 2 median, 3 pca
%  ilndocn: 0 ('o'cn), 1 ('l'nd), 2 (ocn + lnd), 3 (no condition)
  %maxnt = 10;% 7 (lnd) 10 (ocn & ocn+lnd)
  variancethr = 0.60; % threshold of total variance ratio to truncate the modes
  subrows = 4; % subplot rows
  subcols = 5; % ceil(maxnt/subrows)
  clustthr = 2;%18; %(nt+2)*parm.nzi
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
  vname = {'spdt','spdq','precc','precl','freq_precc','freq_precl','divmax','LI'};
  display(vname{ivv})
  vdim = [1 1 0 0 0 0 1 0]; % 0 (0D) 1 (1D)
  ylimrange = {[-0.3 18],[-0.3 18],[0 5],[0 5],[0 0.5],[0 .5]};
  lndocnstr = {'_ocn','_lnd','_lndocn','nocond'};
  ilndocn = ilndocn + 1; % +1 for actual index
  stats = {'mean','median','pca'};
  eval(sprintf('load(''mcs_stats%s.mat'',''%s'',''nlndocnCl'',''ntvec2'');',lndocnstr{ilndocn},vname{ivv}));
  eval(sprintf('vtmp = %s;',vname{ivv}));
  load('mcs_cluster_parm')
  if vdim(ivv) >= 1
    cmapname = {'cmap_dt1','cmapdt1','cmap(1)','','','','cmap(1)'};
    caxislim = {[-5 20]/86400,[-5 20],[-0.005 0.005],[],[],[],[-5e-4 5e-4]};
    eval(sprintf('load(''cmaps'',''%s'');',cmapname{ivv}))
    eval(sprintf('cmapmat = %s;',cmapname{ivv}))
  else
  end
  for is = 1:numel(nlndocnCl)-2
    nt=ntvec2(is);
    if nlndocnCl{nt} >= clustthr % if clusters exists
      timerange = [0:3:3*(nt-1)];
      if ismember(istat,[1 2])
        eval(sprintf('vtmp2{nt} = vtmp.%s{nt};',stats{istat}));
      elseif istat == 3
        varitmp = vtmp.variance{nt};
        pc_meantmp = vtmp.pc_mean{nt};
        svectmp = vtmp.svec{nt};
        nm=1;
        varitmp1 = sum(varitmp(1)/sum(varitmp));
        while varitmp1 < variancethr % set truncation modes according to total variance ratio 
          nm = nm + 1;
          varitmp1 = sum(varitmp(1:nm)/sum(varitmp));
        end
        truncm(nt) = nm; 
        sumvar = 0; 
        for im = 1:truncm(nt); 
          sumvar = sumvar + pc_meantmp(im)*svectmp(:,im); % summing the truncated variable
        end 
        vtmp2{nt}=sumvar; 
      end
      if ~exist('ha')  
%        ha(nt-1) = subplot(subrows,subcols,nt-1);
        if nt==3
          [ha, pos] = tight_subplot(subrows, subcols, [0.08 0.01], 0.05, 0.1)%marg_h, marg_w);
        end
      else
%        subplot(subrows,subcols,nt-1,ha(nt-1)); hold on;
      end
      switch vdim(ivv)
      case 0
        axes(ha(is))
        yyaxis right
        %plot(vtmp{nt},'r');
        if ismember(ivv,[3 5])
          %bar(vtmp2{nt},0.4,'facecolor','none','edgecolor',[.25 .25 .25]);
          plot(timerange,vtmp2{nt},'r')
        elseif ismember(ivv,[4 6])
          plot(timerange,vtmp2{nt},'b');hold on
          plot(timerange,vtmp1{nt},'r')
%          hb = bar([vtmp2{nt}',vtmp1{nt}'],'stacked','barwidth',.4); %,'edgecolor',[.25 .25 .25]);
%          set(hb(1),'facecolor','b')
%          set(hb(2),'facecolor','r')
        end
        set(ha(is),'XTick',timerange)
        ylim(ylimrange{ivv});
        %set(ha(is),'YTick',linspace(ylimrange{ivv}(1),ylimrange{ivv}(2),6),'ycolor','r')
        set(ha(is),'YTick',[.25 .5],'ycolor',[.25 .25 .25])
        %set(ha(is),'YTicklabel',[.5 1],'ycolor',[.25 .25 .25])
        if ~ismember(is,subcols:subcols:subrows*subcols)
          set(ha(is),'YTickLabel',[])
        end
        %set(ha(is),'xtick',[1:nt])
        %xlim([1 nt])
        %set(ha(is),'xtickLabel',3*[0:nt])
      case 1
        axes(ha(is))
        %yyaxis left
        contourf(timerange,parm.zint/1000,reshape(vtmp2{nt},parm.nzi,nt),80,'linestyle','none');
        set(ha(is),'YTickLabel','')
        set(ha(is),'YTick','')
        if ismember(is,1:subcols:subrows*subcols) % if not the first col subplot then remove the ticks
          set(ha(is),'YTick',round([parm.zint(1:2:end)/1000],1))
%          set(ha(is),'YTickLabel',round([parm.zint(1:2:end)/1000],1),'ycolor','k')
        end
        ylim(ylimrange{ivv})
        caxis(caxislim{ivv}); 
%        if mod(nt-1,subcols)==0 % only use colorbar at the right edges cells
%          colorbar; 
          colormap(cmapmat);
%        end
        set(ha(is),'xtick',timerange)
        if istat == 3
          title(['n = ' num2str(nlndocnCl{nt}) '   m = ' num2str(truncm(nt)) '   var% = ' num2str(round(varitmp1,2)*100) '%'])
        else 
          title(['n = ' num2str(nlndocnCl{nt}) ])
        end
        %yyaxis right
      end
    end
  end
  if ismember(ivv,[3 5])
    vtmp1 = vtmp2;
  end
  set(gcf,'color','w')
