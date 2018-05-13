function [ivv] = mcs_cluster_stats(diro,ilndocn,season)
%function mcs_cluster_stats(vname,vdim,camcrm{icc})
  % Purpose: Calc mean, median or eof for cluster samples in 1Dspace-time (heating) or 0Dspace-time (precip). 
  %   The input data random variable dimension is the space and time index
  %   (e.g., z-t), and the sample dimension is the clusters. If the space index
  %   is just a point (e.g., surface rain rate), then the r.v. dimension is
  %   just time.
  % Input Arguments:
  %   vname: 'spdt', 'dtcond', 'spdq', 'prect',  etc 
  %   vdim: 0, 1
  %   camcrm: 'CAM', 'CRM'
  %    
% mcs_cluster_stats(1,0)
%ivv = 2 % 2 spdt 3 dtcond 4 spdq 5 prect
%ivv = 1 % precc precl divmax_crm LI_crm FFTke_crm
%ivv = 2 % freq_precc freq_precl du03km dv03km 
%ilndocn = 2 % 0 ('o'cn), 1 ('l'nd), 2 (ocn + lnd), 3 (no condition)
  doshortclust = 1; % cluster lifetime of 3 hrs
  run ~/scripts/matlab/startup.m
  minnt = 3
  maxnt_intrp = 20;
  maxnt = 20; % aliasing issue if maxnt_intrp< maxnt
  %maxnt = 22
%  vnamei={'dU03km_stormrel'};
%  vnameo={'dU03km_stormrel'};;
%  vdim = zeros(1,1);
%  vnamei = {'B'};
%  vnameo = {'B'};
%  vdim = [zeros(1,1)];
  vnamei = {'ddq_cam_sp'};
  vnameo = {'ddq_cam_sp'};
  vdim = [ones(1,1)];
%  vnamei = {'ddt_cam_sp','ddq_cam_sp'};
%  vnameo = {'ddt_cam_sp','ddq_cam_sp'};
%  vdim = [ones(1,2)];
  convstrat=[];
%  vnamei={'precc','precl','frac_precc','frac_precl'};
%  vnameo={'precc','precl','frac_precc','frac_precl'};
%  vdim = [zeros(1,4)];
%  vnamei={'LIave','LImax','LImin','Bave','Bmax','du03km_cam','dv03km_cam','dU03km_cam','du03km'};
%  vnameo={'LIave','LImax','LImin','Bave','Bmax','du03km_cam','dv03km_cam','dU03km_cam','du03km'};
%  vdim = [zeros(1,9)];
%  vnamei={'LIave','LImax','Bave','Bmax','dU03km_cam'};
%  vnameo={'LIave','LImax','Bave','Bmax','dU03km_cam'};
%  vdim = [zeros(1,5)];
%  vnamei={'du03km_cam','du03km_cam_abs'};
%  vnameo={'du03km_cam','du03km_cam_abs'};;
% 'ddt_cam_sp','ddq_cam_sp',
%vnamei={'camdq','mpdq','macpdq'};
%vnameo={'camdq','mpdq','macpdq'};
%vdim = [ones(1,3)];
%{
vnamei={'camdt','camdq','spdt','spdq','zmdt','evaptzm','cmfdt','mpdt','macpdt','zmdq','evapqzm','cmfdq','mpdq','macpdq',...
        'w_crm','w_crm','qv_crm','qv_crm','qc_crm','qc_crm','qi_crm','qi_crm','qr_crm','qr_crm','frac_precc','frac_precl',...
  'precc','precl'};
vnameo={'camdt','camdq','spdt','spdq','zmdt','evaptzm','cmfdt','mpdt','macpdt','zmdq','evapqzm','cmfdq','mpdq','macpdq',...
        'w_strat','w_conv','qv_strat','qv_conv','qc_strat','qc_conv','qi_strat','qi_conv','qr_strat','qr_conv','frac_precc','frac_precl','precc','precl'};
vdim = [ones(1,14),2*ones(1,10),zeros(1,4)];
%}
%vnamei={'camdt','camdq','spdt','spdq','zmdt','evaptzm','cmfdt','mpdt','macpdt','zmdq','evapqzm','cmfdq','mpdq','macpdq'};
%vnameo={'camdt','camdq','spdt','spdq','zmdt','evaptzm','cmfdt','mpdt','macpdt','zmdq','evapqzm','cmfdq','mpdq','macpdq'};
  convstrat = containers.Map({'w_strat','w_conv','qv_strat','qv_conv','qi_strat','qi_conv','qc_strat','qc_conv','qr_strat','qr_conv','rh_strat','rh_conv'},[0,1,0,1,0,1,0,1,0,1,0,1]); % convstrat('w_strat')=0
%{
  vnamei = {'w_crm','w_crm','qv_crm','qv_crm','qc_crm','qc_crm','qi_crm','qi_crm','qr_crm','qr_crm','rh_crm','rh_crm','u_crm','T_crm','rho_crm','buoy', ... %16
   'spdt','spdq','spmc','spmcup','thetae','relhum','camdt','camdq','ddt_cam_sp','ddq_cam_sp'... %10
   'prect','precc','precl','precc_cam','freq_precc','freq_precl','frac_precc','frac_precl','frac_prec0','frac_prec3','frac_prec5','du03km','dv03km','du03km_abs','dv03km_abs','du03km_cam_abs','flnt','fsnt','B'}; %19
  vnameo = {'w_strat','w_conv','qv_strat','qv_conv','qc_strat','qc_conv','qi_strat','qi_conv','qr_strat','qr_conv','rh_strat','rh_conv','u_crm','T_crm','rho_crm','buoy', ... %16
   'spdt','spdq','spmc','spmcup','thetae','relhum','camdt','camdq','ddt_cam_sp','ddq_cam_sp'... %10
   'prect','precc','precl','precc_cam','freq_precc','freq_precl','frac_precc','frac_precl','frac_prec0','frac_prec3','frac_prec5','du03km','dv03km','du03km_abs','dv03km_abs','du03km_cam_abs','flnt','fsnt','B'}; %19
  vdim = [2*ones(1,16), ones(1,10), zeros(1,19)];
%}
  nvv = numel(vnamei);
  lndocnstr = {'_ocn','_lnd','_lndocn','nocond'}; % removed coast from trackMCS.m!!!!
  ilndocn = ilndocn + 1; % +1 for actual index
%  spcase      = 'F_2000_SPCAM_m2005_3hrly1';
%  spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
%  spcase      = 'F_2000_SPCAM_m2005_3hrly_f09f09_branch';
%  spArchive   = ['/glade2/scratch2/ginochen/archive/matlab/'];
%  diro  = [spArchive spcase '/atm/hist/' season];
  load([diro '/mcs_clusters_1.mat'],'t','t4Cl','nt4Cl','mcsnll','nCl4nt','nCl','nt4ClRowInd','mcsillt4Cl')
  load([diro '/mcs_cluster_parm.mat']);
  load([diro '/mcs_cluster_var/mcs_cluster_ctype']);
  load([diro '/mcs_cluster_var/mcs_cluster_prate']);
  load([diro '/mcs_cluster_var/mcs_cluster_mcsllx4Cl']);
  nz = parm.nz; % 10km=12; 15km=17
  numt = numel(t);
  if ~exist('maxnt','var')
    maxnt = max(nt4Cl)+2;
  end
  if doshortclust
    ntvec2 = [3:maxnt]; % include shortest lifetime (i.e., nt=3) cluster
  else
    ntvec2 = [4:maxnt];
  end
  if strcmp(lndocnstr{ilndocn},'nocond')
    island = zeros(nCl,1);
    disp('setting island to a zero vector, and make sure to set lndocn = 0 to do all clusters without lnd ocn condition')
  else
    load([diro '/mcs_clusters_1.mat'],'island');
  end
  fo = [diro '/mcs_cluster_stats/mcs_stats' lndocnstr{ilndocn} '.mat']; 
  if ~exist(fo,'file')
    disp([ fo ' does not exists, creating a new file...']);
    save(fo,'ntvec2'); 
  else
    disp([ fo ' exists, make sure you want to add more variable statistics, continue running now...']);
    tmp = whos('-file',fo);
    v_exist = {tmp.name};
    ii = 1;
    for iv = 1:nvv
%      if ~ismember(vnameo{iv},v_exist) % only do the stats for variable not in the file
        tmp1{ii} = vnameo{iv}; 
        tmp2(ii) = vdim(iv);
        ii=ii+1;
%      end
    end 
    vnameo = tmp1;
    disp(['Missing variables are: ' vnameo])
    vdim = tmp2;
  end
  ntvec0 = ntvec2-2; % original minimum active mcs lifetime is 1, adding one to both ends become a minimum of 3
  [iclndocn nlndocnCl] = selectlndocnCluster(ntvec0,nt4Cl,nCl4nt,nt4ClRowInd,ilndocn,island,t4Cl,numt); % season head-tail clusters discarded here
  iclndocn2 = cat(2,iclndocn{:}); % cluster index ordered as var.data
  if ~ismember('iclndocn',who('-file',fo)) % if 'iclndocn is not in fo then save it
    save(fo,'iclndocn','iclndocn2','nlndocnCl','-append')
  end
%  if vdim(1)==2 % get the conv and strat zonal indices for just the 2D variable 
    [ixconv ixstrat] = convstratidx(ntvec2,nlndocnCl,iclndocn,ctype,mcsllx4Cl);
%  end
  for ivv = 1:nvv
    disp(vnameo{ivv})
    fi = [diro '/mcs_cluster_var/mcs_cluster_' vnamei{ivv} '.mat'];
    vtmp = loadSingleVar(fi,vnamei{ivv});
    ncc = 0; % cluster counter
    vtmp2 = [];
    eval(sprintf('%s = [];',vnameo{ivv}));
    for it = 1:numel(ntvec2)
      nt2 = ntvec2(it);  
    % nt2 = ntvec2(nCl4nt>1);
tic
      if nlndocnCl{nt2} == 0; continue; end % cannot calc stats with one sample
      [vtmp1] = clusterave(vtmp,vdim(ivv),convstrat,vnameo{ivv},nt2,nlndocnCl,iclndocn,ctype,ixstrat,ixconv);
%      eval(sprintf('%s = mystats(%s,vtmp1,nt2,nlndocnCl);',vname{ivv},vname{ivv})) 
      display('Start lifetime composite statistics calculation, interpolate lifetime into 100%')
      [vtmp2] = mcslifetimeinterp(vtmp2,nt2,maxnt_intrp,vtmp1,vdim(ivv),nz);
toc
    end
    eval(sprintf('%s = mcscomposite_stats(%s,vtmp2,vdim(ivv))',vnameo{ivv},vnameo{ivv}))
    try
      eval(sprintf('save(''%s'',''%s'',''-append'')',fo,vnameo{ivv}))
    catch
      disp('Removing and resaving a variable from the old file')
      eval(sprintf('vars = rmfield(load(''%s''),''%s'');',fo,vnameo{ivv}));
      save(fo,'-struct','vars');
    end
  end
  if vdim(1)==2
    try 
      eval(sprintf('save(''%s'',''%s'',''%s'',''-append'')',fo,'ixstrat','ixconv'))
    catch
      disp('skipping the saving of ''ixstrat'' and ''ixconv''');
    end
  end
  disp('continue run plot_mcs_clust_stats.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ixconv ixstrat] = convstratidx(ntvec2,nlndocnCl,iclndocn,ctype,mcsllx4Cl)
  % Purpose: obtain the convective and stratiform zonal indices for each LS point in each cluster
    for it = 1:numel(ntvec2)
      nt2 = ntvec2(it);
      if nlndocnCl{nt2} <= 1; continue; end % cannot calc stats with one sample
      for iic = 1:nlndocnCl{nt2}
        ic = iclndocn{nt2}(iic);
        for iit = 1:nt2
          for ill=1:size(ctype{ic}{iit},3)
            ixconv{ic}{iit}{ill}  = intersect(find(ctype{ic}{iit}(:,ill)==1),mcsllx4Cl{ic}{iit}{ill}); % conv or strat indices in each mcs
            ixstrat{ic}{iit}{ill} = intersect(find(ctype{ic}{iit}(:,ill)==0),mcsllx4Cl{ic}{iit}{ill}); % conv or strat indices in each mcs
          end
        end
      end
    end


function [vtmp1] = clusterave(vtmp,vdim,convstrat,vnameo,nt2,nlndocnCl,iclndocn,ctype,ixstrat,ixconv)
  % Purpose: Obtain the cluster large-scale spatial ave of a variable at a fixed time
      vtmp1 = []; 
      for iic = 1:nlndocnCl{nt2}
        ic = iclndocn{nt2}(iic);
        for iit = 1:nt2
          switch vdim
          case 0
            vtmp1(iit,iic)     = mean(vtmp{ic}{iit}); % ndims goes to the last index
          case 1
            vtmp1(:,iit,iic)   = mean(vtmp{ic}{iit},2);
          case 2
            %%%%% CONV-STRAT LS average for 2D crm fields %%%%%%%%%%%%%%%%%%%
            vtmp3 = [];
            for ill=1:size(ctype{ic}{iit},3)
%              ixcs = intersect(find(ctype{ic}{iit}(:,ill)==convstrat(vnameo{ivv})),mcsllx4Cl{ic}{iit}{ill}); % conv or strat indices in each mcs
              try  
                if convstrat(vnameo)==0
                    vtmp3(:,ill)   = mean(vtmp{ic}{iit}(ixstrat{ic}{iit}{ill},:,ill),1);
                elseif convstrat(vnameo)==1
                  vtmp3(:,ill)   = mean(vtmp{ic}{iit}(ixconv{ic}{iit}{ill},:,ill),1);
                end 
              catch % if variable doesn't require conditioning on CONV or STRAT then do all zonal average
                vtmp3(:,ill)   = mean(vtmp{ic}{iit}(:,:,ill),1); 
              end
            end
            vtmp1(:,iit,iic)   = mean(vtmp3,2);
            %vtmp1(:,:,iit,iic) = vtmp{ic}{iit}(ix,iz,ill);
%            if iit==1 % first and last index is added when doing mcs_cluster_var.m and mcs_cluster_glueVar.m
%              illc = illtcentroid{ic}(1);
%            elseif iit == nt2
%              illc = illtcentroid{ic}(end);
%            else
%              illc = illtcentroid{ic}(iit-1);
%            end
%            vtmp1(:,:,iit,iic) = vtmp{ic}{iit}(:,:,illc);
          end
        end
      end 

function [vtmp2] = mcslifetimeinterp(vtmp2,nt2,maxnt_intrp,vtmp1,vdim,nz)
  % Purpose: Interpolate the mcs into a 100% lifetime
      to = linspace(0,100,nt2)'; % original time, e.g., for nt2=3, to = [0 50 100], and tq = [ 0    5.2632   10.5263... 94.7368  100.0000]
      tq = linspace(0,100,maxnt_intrp)'; % interpolated time
      if vdim==0
        vtmp2{nt2} = interp1(to,vtmp1,tq,'linear'); % calc percentage lifecycle
      elseif ismember(vdim,[1,2])%ismember(vdim,[1,2])
%        for iz = 1:parm.nzi
        for iz = 1:nz
          vtmp2{nt2}(iz,:,:) = interp1(to,squeeze(vtmp1(iz,:,:)),tq,'linear'); % calc percentage lifecycle
        end
%      elseif vdim==2
%        for ix = 1:parm.nx
%          for iz = 1:iz10km
%            vtmp2{nt2}(ix,iz,:,:) = interp1(to,squeeze(vtmp1(ix,iz,:,:)),tq,'linear'); % calc percentage lifecycle
%          end
%        end
      else
        error('dimension non-existence');
      end


function myv = mcscomposite_stats(myv,vtmp2,vdim)
    id = ndims(vtmp2{3}); % start from nt2=3, id = sample dimension
    vtmp = cat(id,vtmp2{:}); % concat cell into array according to id-dim
    nc = size(vtmp,id);
    myv.mcscomposite.data = vtmp;
    myv.mcscomposite.mean = nanmean(vtmp,id); 
    myv.mcscomposite.median = nanmedian(vtmp,id); 
    myv.mcscomposite.std = nanstd(vtmp,[],id); 
    doecon = 1;
    dosignflip=0
    [myv.mcscomposite.svec myv.mcscomposite.pcvariance myv.mcscomposite.pc myv.mcscomposite.pc_mean] = svdVar(reshape(vtmp,numel(vtmp)/nc,nc)',doecon,dosignflip); % varin(nsamp,nvar)


function [iclndocn nlndocnCl] = selectlndocnCluster(ntvec0,nt4Cl,nCl4nt,nt4ClRowInd,ilndocn,island,t4Cl,numt)
  %Purpose: obtain cluster indices assoc to lnd/ocn condition
  %Input: 
  %  ntvec0: original active MCS lifetime, without the growth and decay extra two time indices
  %Output: 
  %  iclndocn: cluster indices assoc to lnd/ocn condition
  %  nlndocnCl: total number of clusters for cluster lifetime nt 
  ii = 1; % discarded cluster counter
  for nt = ntvec0
      iiic = 1;
      for iic = 1:nCl4nt(nt)
        ic = nt4ClRowInd{nt}(iic); % ic-cluster with lifetime nt = i
        %if any(ismember(t4Cl(ic,:),[1 numt])) | nt4Cl(ic)==1; continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
        if any(ismember(t4Cl(ic,:),[1 numt])); ii=ii+1; continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step; the removed clusters are iCends
        if island(ic)==ilndocn-1  % qualify for ocean index or default no land ocean condition
          iclndocn{nt+2}(iiic) = ic;
          iiic = iiic+1; % counter for lnd ocn conditioned cluster
        end
      end   
      nlndocnCl{nt+2} = iiic-1; % total count of lnd ocn conditioned cluster 
  end
  disp(['total discarded clusters = ' num2str(ii)]);

function vtmp = loadSingleVar(fi,vname);
% load vname into vtmp without knowing the name of vname{ivv}
    vtmp = load(fi,vname);
    foo = fieldnames(vtmp);
    vtmp = vtmp.(foo{1});

function [myv] = mystats(myv,vtmp1,nt,nlndocnCl)
  % myv: my variable 
  disp('begin median calc')
  myv.median{nt} = median(vtmp1,ndims(vtmp1));
  myv.mean{nt} = mean(vtmp1,ndims(vtmp1));
  disp('end median calc')
  disp(['begin pca, mean, variance for time length ' num2str(nt) ])
  [myv.svec{nt} myv.pcvariance{nt} myv.pc{nt} myv.pc_mean{nt}] = svdVar(reshape(vtmp1,numel(vtmp1)/nlndocnCl{nt},nlndocnCl{nt})'); % varin(nsamp,nvar)
  disp(['finish pca for time length ' num2str(nt) ])
