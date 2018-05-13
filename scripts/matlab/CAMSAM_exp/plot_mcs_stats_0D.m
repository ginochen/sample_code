%  figure; 
  istat = 2;
  imm = 2;
  icc = 1;
  ylimrange = [0 20];
  vname = 'prect' % prect:prect, dcq:spdq for CAM and CRM are the same 
  stats = {'meanmedian','pca'};
  meanmedian = {'mean','median'};
  camcrm = {'CAM','CRM'};
  eval(sprintf('load mcs_stats_%s_%s_%s.mat;',stats{istat},camcrm{icc},vname));
  subrows = 3; % subplot rows
  minnt = 2;
  maxnt = 12;
  truncm = 1;
  for nt=minnt:maxnt;
      switch istat
      case 1
        eval(sprintf('vtmp{nt} = %s.%s{nt};',vname,meanmedian{imm}));
      case 2 % truncated svd sum
        sumvar = 0; 
        for im = 1:truncm; 
          eval(sprintf('sumvar = sumvar + (%s.svec{nt}(:,im)*%s.pc_mean{nt}(im));',vname,vname)); % summing the truncated variable
        end; 
        vtmp{nt}=sumvar; 
      end
    if ~isempty(vtmp{nt}) % if mean not empty then plot, otherwise skip
      subplot(subrows,ceil(maxnt/subrows),nt,ha(nt));
      yyaxis right
      plot(vtmp{nt},'k');
      ylim(ylimrange);
      xlim([1 nt])
    end
  end
