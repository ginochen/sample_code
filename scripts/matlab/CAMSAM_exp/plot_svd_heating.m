% heating eof & mass flux
  load mcs_clusters; % pick out lengthly mcs clusters to plot
  dtcam1D = [];
  dtcrm1D = [];
% create index variables 
% assemble var into var(randvar,sample)
  fname = 'F_2000_SPCAM_m2005_3hrly1_MCS';
  ii = 1;
  mday = [31 28 31 30 31 30 31 31 30 31 30 31]; % # of days in a month
  % check if .mat file exists
  for iy = 1:1:1
    for im = 2:1:2 
      for id = 2:1:mday(im)
        for is = 0:10800:75600
          time=num2str(sprintf('%04g-%02g-%02g-%05g',iy,im,id,is));
          if exist([fname '.' time '.mat'])
            t{ii} = time;
            ii=ii+1;
          end
        end
      end
    end
  end
  if ( ~exist('t') )
    error('no file to track')
  end
  load([fname '.' t{1} '.mat'],'parm','idv')
  for it = 1:numel(t)
    eval(sprintf('load(''F_2000_SPCAM_m2005_3hrly1_MCS.%s.mat'',''varCAM1D'');',t{it}));
    dtcam1D = cat(2,dtcam1D,varCAM1D(1:parm.nzi,:,idv.varCAM1D.dtcond)); % the index of columns is the same as index of rows of mcslltC
    dtcrm1D = cat(2,dtcrm1D,varCAM1D(1:parm.nzi,:,idv.varCAM1D.spdt));
  end
  
  N=3 % 2 steps is only 3hours apart, not 6 hours
  ir = find(nt4Cl>=N); % cluster with # of time steps greater >= N
% obtain singular vector 
  %[svec_cam sval_cam pc_cam] = svdVar(dtcam1D(:,isnan(dtcam1D(1,:))~=1)); % var(randvar,sample)
  %[svec_crm sval_crm pc_crm] = svdVar(dtcrm1D(:,isnan(dtcrm1D(1,:))~=1)); % var(randvar,sample)
  [svec_cam sval_cam pc_cam] = svdVar(dtcam1D(:,ir)); % var(randvar,sample)
  [svec_crm sval_crm pc_crm] = svdVar(dtcrm1D(:,ir)); % var(randvar,sample)
  corr(pc_cam(:,1),pc_crm(:,1))
  sval_ratio_cam = sval_cam/sum(sval_cam)
  sval_ratio_crm = sval_crm/sum(sval_crm)
% subplot 1st 3 eigvec in vert plane
  neof=3;
  zint=parm.zint/1000; % m to km
% mean plot
  ylim_vec=[-5 22; -2 3; -1 1.5];
  figure
  for ieof=1:neof
    subplot(1,neof,ieof)
    plot(zint,svec_crm(:,ieof)*mean((pc_crm(:,ieof))) ,'color',[.5 .5 .5]);hold on;
    plot(zint,svec_crm(:,ieof)*median((pc_crm(:,ieof))),'--','color',[.5 .5 .5] );hold on;
    plot(zint,svec_cam(:,ieof)*mean((pc_cam(:,ieof))), 'r');
    plot(zint,svec_cam(:,ieof)*median((pc_cam(:,ieof))), 'r--');
    ylim(ylim_vec(ieof,:))
    camroll(90);set(gca,'Xdir','reverse')
  end
  
  set(gcf,'color','w')

