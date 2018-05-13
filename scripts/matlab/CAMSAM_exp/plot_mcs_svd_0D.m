  vname = 'prect';
  eval(sprintf('load mcs_stats_pca_%s.mat;',vname));
  subrows = 3; % subplot rows
  minnt = 2;
  maxnt = 10;
  truncm = 1; % truncated at truncm-mode
%  caxislim = [-2 20];
%  cmapmat = cmap_spdt_nt9;
cmapmat = cmap(1);
  figure
  for nt=minnt:maxnt; 
    sumvar = 0; 
    for im = 1:truncm; 
      eval(sprintf('sumvar = sumvar + (%s.svec{nt}(:,im)*%s.pc_mean{nt}(im));',vname,vname)); % summing the truncated variable
    end; 
    vtmp{nt}=sumvar; 
    subplot(subrows,ceil(maxnt/subrows),nt);
    plot(1:nt,vtmp{nt});
    caxis(caxislim); 
    colorbar; 
    colormap(cmapmat);
  end
