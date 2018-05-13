% mass flux
timesteps = 41:10:401;
ucam1D_cmfmc = [];
ucam1D_cmfmcdzm = [];
ucrm1D = [];

% assemble var into var(randvar,sample)
for it = timesteps
  load(['var_PC1_' num2str(it) '_mcc.mat'],'varCAM1D','parm')
  % create index variables 
  if it==timesteps(1)
    idv = varindex({'varCAM1D'},parm)
  end
  ucam1D_cmfmc = cat(2,ucam1D_cmfmc,varCAM1D(1:parm.nzi,:,idv.varCAM1D.cmfmc));
  ucam1D_cmfmcdzm = cat(2,ucam1D_cmfmcdzm,varCAM1D(1:parm.nzi,:,idv.varCAM1D.cmfmcdzm));
  ucrm1D = cat(2,ucrm1D,varCAM1D(1:parm.nzi,:,idv.varCAM1D.spmcup));
end
ucam1D = ucam1D_cmfmc + ucam1D_cmfmcdzm;
% obtain singular vector 
[svec_cam sval_cam pc_cam] = svdVar(ucam1D(:,isnan(ucam1D(1,:))~=1)); % var(randvar,sample)
[svec_crm sval_crm pc_crm] = svdVar(ucrm1D(:,isnan(ucrm1D(1,:))~=1)); % var(randvar,sample)
corr(pc_cam(:,1),pc_crm(:,1))
sval_ratio_cam = sval_cam/sum(sval_cam)
sval_ratio_crm = sval_crm/sum(sval_crm)
% subplot 1st 3 eigvec in vert plane
neof=3;
zint=parm.zint/1000; % m to km
% median plot
figure
for ieof=1:neof
  subplot(1,neof,ieof)
  plot(zint,svec_crm(:,ieof)*mean((pc_crm(:,ieof))), 'color', [0.5 0.5 0.5] );hold on;
  plot(zint,svec_crm(:,ieof)*median((pc_crm(:,ieof))),'linestyle','--', 'color', [0.5, 0.5, 0.5] );
  plot(zint,svec_cam(:,ieof)*mean((pc_cam(:,ieof))), 'r' );
  plot(zint,svec_cam(:,ieof)*median((pc_cam(:,ieof))),'r--' );camroll(90);set(gca,'Xdir','reverse')
end



function idv = varindex(str,parm)
  % idv: object variable for variable index
  % example: 
  %   % the index for SST variable
  %   idx.SST = 1
  %   % the str is thus 'SST'
  %   % parm saves the variable in index order
  %   % so suppose two variables 'SST' and 'SLP'
  %   parm = {'SST', 'SLP'}
  for ii=1:numel(str)
    evalc(sprintf('nv = numel(parm.%s)',str{ii}));
    evalc(sprintf('strvar = parm.%s',str{ii}));
    for iv = 1:nv
      evalc(sprintf('idv.%s.%s = %d',str{ii},strvar{iv},iv));
    end
  end
