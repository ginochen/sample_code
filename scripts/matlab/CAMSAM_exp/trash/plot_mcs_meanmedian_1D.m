  figure; 
  imm = 1;
  icc = 2;
  vname = 'spdt' % prect:prect, dcq:spdq for CAM and CRM are the same 
  meanmedian = {'mean','median'};
  camcrm = {'CAM','CRM'};
  eval(sprintf('load mcs_stats_meanmedian_%s_%s.mat;',camcrm{icc},vname));
  subrows = 3; % subplot rows
  minnt = 2;
  maxnt = 22;
  truncm = 1;
  for nt=minnt:maxnt;
    if eval(sprintf('~isempty(%s.%s{nt});',vname,meanmedian{imm})) % if mean not empty then plot, otherwise skip
      subplot(subrows,ceil(maxnt/subrows),nt);
      eval(sprintf('contour(1:nt,parm.zint,reshape(%s.%s{nt},parm.nzi,nt),30)',vname,meanmedian{imm}));
      xlim([1 nt])
    end
  end
