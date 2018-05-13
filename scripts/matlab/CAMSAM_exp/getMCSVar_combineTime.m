  cd /projects/rsmas/kirtman/gchen/cesm_spcam/archive/F_2000_SPCAM_m2005_3hrly1/atm/hist
  fname = 'F_2000_SPCAM_m2005_3hrly1_MCSVar'; 
  load('mcs_clusters','t','t4Cl','nt4Cl','mcsnll','nCl4nt','nt4ClRowInd','mcsillt4Cl');
  eval(sprintf('load(''%s.%s.mat'',''parm'',''idv'');',fname,t{1}));
  for it = 1:numel(t)
    tic
    eval(sprintf('load(''%s.%s.mat'');',fname,t{it}))
    varCAM0D_all{it} = varCAM0D; 
    varCRM0D_all{it} = varCRM0D; 
    varCAM1D_all{it} = varCAM1D; 
    varCRM1D_all{it} = varCRM1D; 
    varCRM2D_all{it} = varCRM2D; 
    stabilityIndex_all{it} = stabilityIndex;
    divmax_crm_all{it} = divmax_crm;
    enstrophy_crm_all{it} = enstrophy_crm;
    FFTke_crm_all{it} = FFTke_crm;
    mcs_all{it} = mcs;
    clear varCAM0D varCRM0D varCAM1D varCRM1D varCRM2D stabilityIndex divmax_crm enstrophy_crm FFTke_crm mcs
    toc
  end
  display('begin saving var all time')
  save('F_2000_SPCAM_m2005_3hrly1_MCSVarMAM','varCAM0D_all','varCRM0D_all','varCAM1D_all','varCRM1D_all','stabilityIndex_all','divmax_crm_all','enstrophy_crm_all','FFTke_crm_all','mcs_all');
  save('F_2000_SPCAM_m2005_3hrly1_MCSVar2DMAM','varCAM2D_all','-v7.3');
  display('finish saving var all time')
