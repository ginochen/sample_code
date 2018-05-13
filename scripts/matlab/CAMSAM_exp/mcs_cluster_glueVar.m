function [vf] = mcs_cluster_glueVar(fmcs,fparm,diro)
% Purpose: combine the data created from mcs_cluster_ext.m
%   the for/if loops are exactly the same as mcs_cluster_ext.m
%close all force % close all mat/nc file pointer opened
donew = 0 % 1 glue new variable, 0 glue all variables
if ~exist('donew')
  donew = 0 % do all variables
end
%spcase      = 'F_2000_SPCAM_m2005_3hrly1';
%spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
%spcase      = 'F_2000_SPCAM_m2005_3hrly_f09f09_branch';
%spArchive   = ['/glade2/scratch2/ginochen/archive/matlab'];
%spHist.dir  = [spArchive spcase '/atm/hist/' season];
%diro = [spHist.dir '/mcs_cluster_var/'];
load(fmcs); % this file supposed to be linked to current dir
load(fparm)
%load([diro '/mcs_cluster_mcsllx4Cl'])
if donew
%  vnamenew = {'mcsprec','total_precc','total_precl','ictop','Tctop','Tctopmin'};
  vnamenew = {'mcsillt4Cl','frac_prec'};
end
% vnames = {'u_cam','du03km_cam','du03km_cam_abs'};
% vnames = {'B','Bave','Bmax','LI','LIave','LImax','LImin','du03km_cam','dv03km_cam','dU03km_cam','du03km'};
% vnames = {'prect'};
% vnames = {'precc','precl'};
% vnames = {'u_storm','v_storm','u_stormrel','v_stormrel','dU03km_stormrel'};
% vnames = {'u_cami','v_cami'};
% vnames = {'u_cam','v_cam'};
% vnames = {'ilat','ilon'};
%   vnames = {'precc','precl','frac_precl','frac_prec0','frac_prec3','frac_prec5',...
%   'frac_precc','freq_precc','freq_precl'}
%else % specify new variables 
%  vnames = {'buoy','B','du03km','dzm2d','dzi2d','dzi1d','camdt','camdq','cmfdt',...
%   'cmfdq','ctype','evapqzm','evaptzm','flnt','frac_precl','frac_prec0','frac_prec3','frac_prec5',...
%   'frac_precc','freq_precc','freq_precl','fsnt','ictop','LI','mcsillt4Cl','mcsllx4Cl','macpdt','macpdq','mcsprec','mpdt','mpdq','prate','precc','precl',...
%   'phis','rho2d','rho1d','spdt','spdq','qv_crm','qc_crm','qi_crm','qr_crm','T_crm','Tctop','Tctopmin','total_precc','total_precl',...
%   'w_crm','z1d','z2d','z3','zmdt','zmdq'};
%  vnames = {'mpdt','mpdq'};
%  vnames = {'w','pv','theta','mcsillt4Cl'};
%  vnames = {'u_cam','v_cam','theta0','thetaa','thetaz0','rho0',};
%  vnames = {'pv','pv0','pva'};
%end
% determine if data from mcs_cluster_ext are complete
nt = numel(t);
for it=2:nt-1
  if ~mcsn(it); continue; end
  if ~exist([diro 'mcs_cluster_var.' num2str(t{it}) '.mat'],'file')
    error(['file mcs_cluster_var.' num2str(t{it}) '.mat doesn''t exist, exiting...']) 
  end
end

% main combining program
vf=[];
for it=2:nt-1 % approx 756 time index
    if ~mcsn(it); continue; end
tic
    load([diro 'mcs_cluster_var.' t{it} '.mat'],'vo'); disp(t{it})
    if ~exist('dv','var'); [vnames dv] = getvardim(vo,nCl); end %only do this once
%    for iv = 1:numel(dv)
%      eval(sprintf('vtmp = vo.%s;',vnames{iv}));
      for ic = 1:nCl
%        if any(ismember(t4Cl(ic,:),[1 nt])) | nt4Cl(ic)==1; continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
        if any(ismember(t4Cl(ic,:),[1 nt])); continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
        if (it == mcsillt4Cl{ic}{1}(1,3)) % cluster's starting time index
          for ill=1:mcsnll(ic,1) % go through the ilon-ilat at iit time for cluster ic
%            [vf] = combineall(ic,1,ill,iv, vtmp, vo, vf, vnames, dv, parm, donew);
            [vf] = combineall(ic,1,ill,vo, vf, vnames, dv, parm, donew);
            [vf] = combineall(ic,2,ill,vo, vf, vnames, dv, parm, donew);
            if (it == mcsillt4Cl{ic}{end}(1,3)) % for short lived clusters with only one time index 
              [vf] = combineall(ic,3,ill,vo, vf, vnames, dv, parm, donew);
            end
          end
        elseif (it < mcsillt4Cl{ic}{end}(1,3) & it > mcsillt4Cl{ic}{1}(1,3))% it's in between start and end of a cluster
          itt = find(t4Cl(ic,:)==it);
          for ill=1:mcsnll(ic,itt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
            [vf] = combineall(ic,itt+1,ill,vo, vf, vnames, dv, parm, donew);
          end
        elseif (it == mcsillt4Cl{ic}{end}(1,3)) % cluster's ending time index
          ntt = nt4Cl(ic);
          for ill=1:mcsnll(ic,ntt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
            [vf] = combineall(ic,ntt+1,ill,vo, vf, vnames, dv, parm, donew);
            [vf] = combineall(ic,ntt+2,ill,vo, vf, vnames, dv, parm, donew);
          end
        end % if loop of calc variables 
      end % for loop of cluster indices
%    end
toc
    close all force
end % for loop of time indices
if donew
%  nvv = numel(vnames);
%else
%  nvv = numel(vnames)+numel(vnamenew);
%  vnames = [vnames,vnamenew'];
  vnames = vnamenew';
end
nvv = numel(vnames);
for ivv = 1:nvv
  eval(sprintf('%s = vf.%s;',vnames{ivv},vnames{ivv}));
%  eval(sprintf('%s = vf{ivv};',vnames{ivv}));
  eval(sprintf('save(''%s/mcs_cluster_%s.mat'',''%s'',''-v6'');',diro,vnames{ivv},vnames{ivv}));
  disp(['finished saving ' vnames{ivv}])
end
disp('continue to run mcs_cluster_stats(ilndocn,season)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vname dv] = getvardim(vo,nCl)
% dv: dim vector of variables
vname = fieldnames(vo);
for ic=1:nCl
  if numel(vo.spdt{ic})>0 & ~isempty(vo.spdt{ic}{1})
    for iv=1:numel(vname)
      dv(iv) = sum(eval(sprintf('size(vo.%s{ic}{1});',vname{iv}))>1); 
    end
    break
  end
end

function  [vf] = combineall(ic,iitt,ill,vo, vf, vn, dv, parm, donew)
%function  [vf] = combineall(ic,iitt,ill, vo, vf, vn, dv, parm, donew)
%donew = 1;% 1 glue new variable, 0 glue all variables

if donew
  for ix=1:parm.nx
    if any(squeeze(vo.qi_crm{ic}{iitt}(ix,:,ill))>1e-6)
      vf.ictop{ic}{iitt}(ix,ill) = max(find(squeeze(vo.qi_crm{ic}{iitt}(ix,:,ill))>1e-6)); % save cloud-top vertical index for ix point
      vf.Tctop{ic}{iitt}(ix,ill) = vo.T_crm{ic}{iitt}(ix,vf.ictop{ic}{iitt}(ix,ill),ill);
    else
      vf.ictop{ic}{iitt}(ix,ill) = NaN;
      vf.Tctop{ic}{iitt}(ix,ill) = NaN;
    end
  end
  vf.Tctopmin{ic}{iitt}(ill) = min(vf.Tctop{ic}{iitt}(:,ill));
  vf.mcsprec{ic}{iitt}(ill) = nansum(vo.prate{ic}{iitt}(vo.mcsllx4Cl{ic}{iitt}{ill},ill)); 
  vf.total_precc{ic}{iitt}(ill) = vo.frac_precc{ic}{iitt}(ill)*vf.mcsprec{ic}{iitt}(ill);  
  vf.total_precl{ic}{iitt}(ill) = vo.frac_precl{ic}{iitt}(ill)*vf.mcsprec{ic}{iitt}(ill);  
%}
%{
switch dv(iv)
case 0  
  vf{iv}{ic}{iitt}(ill) = vtmp{ic}{iitt}(ill);
case 1
  if ismember(vn{iv},{'frac_prec','mcsillt4Cl'})
    vf{iv}{ic}{iitt}(:,ill) = vtmp{ic}{iitt}(ill,:);
  else
    vf{iv}{ic}{iitt}(:,ill) = vtmp{ic}{iitt}(:,ill);
  end
case 2
  vf{iv}{ic}{iitt}(:,:,ill) = vtmp{ic}{iitt}(:,:,ill);
end
%if ~donew
if iv==numel(dv) & donew 
  vf{numel(dv)+1}{ic}{iitt}(ill) = nansum(vo.prate{ic}{iitt}(vo.mcsllx4Cl{ic}{iitt}{ill},ill));  % mcsprec
  vf{numel(dv)+2}{ic}{iitt}(ill) = vo.frac_precc{ic}{iitt}(ill)*vf{numel(dv)+1}{ic}{iitt}(ill); % total_precc
  vf{numel(dv)+3}{ic}{iitt}(ill) = vo.frac_precl{ic}{iitt}(ill)*vf{numel(dv)+1}{ic}{iitt}(ill); % total_precl
%}
%{
  vf.prate{ic}{iitt}(:,ill)    =     vo.prate{ic}{iitt}(:,ill)    ;
  vf.precc{ic}{iitt}(ill)      =     vo.precc{ic}{iitt}(ill)      ;
  vf.frac_precl{ic}{iitt}(ill) =     vo.frac_precl{ic}{iitt}(ill) ;
  vf.frac_prec0{ic}{iitt}(ill) =     vo.frac_prec{ic}{iitt}(ill,1);
  vf.frac_prec3{ic}{iitt}(ill) =     vo.frac_prec{ic}{iitt}(ill,2);
  vf.frac_prec5{ic}{iitt}(ill) =     vo.frac_prec{ic}{iitt}(ill,3);
  vf.frac_precc{ic}{iitt}(ill) =     vo.frac_precc{ic}{iitt}(ill) ;
  vf.freq_precc{ic}{iitt}(ill) =     vo.freq_precc{ic}{iitt}(ill) ;
  vf.freq_precl{ic}{iitt}(ill) =     vo.freq_precl{ic}{iitt}(ill) ;
  vf.precl{ic}{iitt}(ill)      =     vo.precl{ic}{iitt}(ill)      ;
  vf.u_cam{ic}{iitt}(:,ill)   =     vo.u_cam{ic}{iitt}(:,ill)   ;
  vf.v_cam{ic}{iitt}(:,ill)   =     vo.v_cam{ic}{iitt}(:,ill)   ;
  vf.u_cami{ic}{iitt}(:,ill)   =     vo.u_cami{ic}{iitt}(:,ill)   ;
  vf.v_cami{ic}{iitt}(:,ill)   =     vo.v_cami{ic}{iitt}(:,ill)   ;
  vf.ilat{ic}{iitt}(ill) =  vo.ilat{ic}{iitt};
  if ill==1 % only use this once
    vf.u_storm{ic}{iitt}    = nanmean(nanmean(vo.u_cami{ic}{iitt}(1:parm.i6km,:)));
    vf.v_storm{ic}{iitt}    = nanmean(nanmean(vo.v_cami{ic}{iitt}(1:parm.i6km,:)));
  end
  vf.u_stormrel{ic}{iitt}(:,ill) = vo.u_cami{ic}{iitt}(:,ill)- vf.u_storm{ic}{iitt};
  vf.v_stormrel{ic}{iitt}(:,ill) = vo.v_cami{ic}{iitt}(:,ill)- vf.v_storm{ic}{iitt};
%  vf.dU03km_stormrel{ic}{iitt}(ill) = sqrt(diff(vf.u_stormrel{ic}{iitt}([1 parm.i3km],ill))^2 +  diff(vf.v_stormrel{ic}{iitt}([1 parm.i3km],ill))^2); % this is same as dU03km_cam
  vf.precc{ic}{iitt}(:,ill)    = vo.precc{ic}{iitt}(ill)*parm.nx/nansum(vo.ctype{ic}{iitt}(:,ill)==1);
  vf.precl{ic}{iitt}(:,ill)    = vo.precl{ic}{iitt}(ill)*parm.nx/nansum(vo.ctype{ic}{iitt}(:,ill)==0);
  vf.prect{ic}{iitt}(:,ill)    = sum(vo.prate{ic}{iitt}(:,ill));
   vf.LImin{ic}{iitt}(ill)     =     min(vo.LI{ic}{iitt}(:,ill))     ;
   vf.B{ic}{iitt}(:,ill)        =     vo.B{ic}{iitt}(:,ill)          ;
   vf.Bave{ic}{iitt}(ill)        =     vo.Bave{ic}{iitt}(ill)          ;
   vf.Bmax{ic}{iitt}(ill)        =     vo.Bmax{ic}{iitt}(ill)          ;
   vf.LI{ic}{iitt}(:,ill)        =     vo.LI{ic}{iitt}(:,ill)          ;
   vf.LImax{ic}{iitt}(:,ill)     =     vo.LImax{ic}{iitt}(ill)          ;
   vf.LIave{ic}{iitt}(:,ill)     =     vo.LIave{ic}{iitt}(ill)          ;
   vf.du03km_cam{ic}{iitt}(ill)  =     vo.du03km_cam{ic}{iitt}(ill)          ;
   vf.dv03km_cam{ic}{iitt}(ill)  =     vo.dv03km_cam{ic}{iitt}(ill)          ;
   vf.dU03km_cam{ic}{iitt}(ill)  =     vo.dU03km_cam{ic}{iitt}(ill)          ;
   vf.du03km{ic}{iitt}(ill)      =     vo.du03km_cam{ic}{iitt}(ill)          ;
%}
else
for iv=1:numel(dv)
  switch dv(iv)
  case 0
   evalc(sprintf('vf.%s{ic}{iitt}(ill)=vo.%s{ic}{iitt}(ill);',vn{iv},vn{iv}));
  case 1 
   if ~ismember(vn{iv},{'frac_prec','mcsillt4Cl'})
     evalc(sprintf('vf.%s{ic}{iitt}(:,ill)=vo.%s{ic}{iitt}(:,ill);',vn{iv},vn{iv}));
   else
     evalc(sprintf('vf.%s{ic}{iitt}(ill,:)=vo.%s{ic}{iitt}(ill,:);',vn{iv},vn{iv}));
   end 
  case 2
   evalc(sprintf('vf.%s{ic}{iitt}(:,:,ill)=vo.%s{ic}{iitt}(:,:,ill);',vn{iv},vn{iv}));
  end
end

%{
   vf.buoy{ic}{iitt}(:,:,ill)   =     vo.buoy{ic}{iitt}(:,:,ill)   ;
   vf.du03km{ic}{iitt}(ill)     =     vo.du03km{ic}{iitt}(ill)     ;
   vf.du03km_cam{ic}{iitt}(ill)   =     vo.du03km_cam{ic}{iitt}(ill)   ;
   vf.du03km_cam_abs{ic}{iitt}(ill)   =     vo.du03km_cam_abs{ic}{iitt}(ill)   ;
   vf.dzm2d{ic}{iitt}(:,:,ill)  =     vo.dzm2d{ic}{iitt}(:,:,ill)  ;
   vf.dzi2d{ic}{iitt}(:,:,ill)  =     vo.dzi2d{ic}{iitt}(:,:,ill)  ;
   vf.dzi1d{ic}{iitt}(:,ill)    =     vo.dzi1d{ic}{iitt}(:,ill)    ;
   vf.camdt{ic}{iitt}(:,ill)    =     vo.camdt{ic}{iitt}(:,ill)    ;
   vf.camdq{ic}{iitt}(:,ill)    =     vo.camdq{ic}{iitt}(:,ill)    ;
   vf.cmfdt{ic}{iitt}(:,ill)    =     vo.cmfdt{ic}{iitt}(:,ill)    ;
   vf.cmfdq{ic}{iitt}(:,ill)    =     vo.cmfdq{ic}{iitt}(:,ill)    ;
   vf.ctype{ic}{iitt}(:,ill)    =     vo.ctype{ic}{iitt}(:,ill)    ;
   vf.evapqzm{ic}{iitt}(:,ill)  =     vo.evapqzm{ic}{iitt}(:,ill)  ;
   vf.evaptzm{ic}{iitt}(:,ill)  =     vo.evaptzm{ic}{iitt}(:,ill)  ;
   vf.flnt{ic}{iitt}(ill)       =     vo.flnt{ic}{iitt}(ill)       ;
   vf.fsnt{ic}{iitt}(ill)       =     vo.fsnt{ic}{iitt}(ill)       ;
   vf.LI{ic}{iitt}(ill)         =     vo.LI{ic}{iitt}(ill)         ;
   vf.macpdt{ic}{iitt}(:,ill)   =     vo.macpdt{ic}{iitt}(:,ill)   ;
   vf.macpdq{ic}{iitt}(:,ill)   =     vo.macpdq{ic}{iitt}(:,ill)   ;
   vf.mcsillt4Cl{ic}{iitt}(ill,1:3) = vo.mcsillt4Cl{ic}{iitt}(ill,1:3);
   vf.mcsllx4Cl{ic}{iitt}{ill} =      vo.mcsllx4Cl{ic}{iitt}{ill}  ;
   vf.mpdt{ic}{iitt}(:,ill)     =     vo.mpdt{ic}{iitt}(:,ill)     ;
   vf.mpdq{ic}{iitt}(:,ill)     =     vo.mpdq{ic}{iitt}(:,ill)     ;
   vf.phis{ic}{iitt}(ill)       =     vo.phis{ic}{iitt}(ill)       ;
   vf.rho2d{ic}{iitt}(:,:,ill)  =     vo.rho2d{ic}{iitt}(:,:,ill)  ;
   vf.rho1d{ic}{iitt}(:,ill)    =     vo.rho1d{ic}{iitt}(:,ill)    ;
   vf.spdt{ic}{iitt}(:,ill)     =     vo.spdt{ic}{iitt}(:,ill)     ;
   vf.spdq{ic}{iitt}(:,ill)     =     vo.spdq{ic}{iitt}(:,ill)     ;
   vf.qv_crm{ic}{iitt}(:,:,ill) =     vo.qv_crm{ic}{iitt}(:,:,ill) ;
   vf.qc_crm{ic}{iitt}(:,:,ill) =     vo.qc_crm{ic}{iitt}(:,:,ill) ;
   vf.qi_crm{ic}{iitt}(:,:,ill) =     vo.qi_crm{ic}{iitt}(:,:,ill) ;
   vf.qr_crm{ic}{iitt}(:,:,ill) =     vo.qr_crm{ic}{iitt}(:,:,ill) ;
   vf.T_crm{ic}{iitt}(:,:,ill)  =     vo.T_crm{ic}{iitt}(:,:,ill)  ;
   vf.u_cam{ic}{iitt}(:,ill)   =     vo.u_cam{ic}{iitt}(:,ill)   ;
   vf.w_crm{ic}{iitt}(:,:,ill)  =     vo.w_crm{ic}{iitt}(:,:,ill)  ;
   vf.z1d{ic}{iitt}(:,ill)      =     vo.z1d{ic}{iitt}(:,ill)      ;
   vf.z2d{ic}{iitt}(:,:,ill)    =     vo.z2d{ic}{iitt}(:,:,ill)    ;
   vf.z3{ic}{iitt}(:,ill)       =     vo.z3{ic}{iitt}(:,ill)       ;
   vf.zmdt{ic}{iitt}(:,ill)     =     vo.zmdt{ic}{iitt}(:,ill)     ;
   vf.zmdq{ic}{iitt}(:,ill)     =     vo.zmdq{ic}{iitt}(:,ill)     ;
%}


%          vf.rh_crm    {ic}{iitt}(:,:,ill)  = vo.rh_crm    {ic}{iitt}(:,:,ill);
%          vf.rho_crm   {ic}{iitt}(:,:,ill)  = vo.rho_crm   {ic}{iitt}(:,:,ill);
%          vf.vort      {ic}{iitt}(:,:,ill)  = vo.vort      {ic}{iitt}(:,:,ill);
%          vf.fftke     {ic}{iitt}(:,:,ill)  = vo.fftke     {ic}{iitt}(:,:,ill);
%          vf.div       {ic}{iitt}(:,ill)    = vo.div       {ic}{iitt}(:,ill);  
%          vf.u_cam     {ic}{iitt}(:,ill)    = vo.u_cam     {ic}{iitt}(:,ill); 
%          vf.v_cam     {ic}{iitt}(:,ill)    = vo.v_cam     {ic}{iitt}(:,ill);     
%          vf.spmc      {ic}{iitt}(:,ill)    = vo.spmc      {ic}{iitt}(:,ill);    
%          vf.spmcup    {ic}{iitt}(:,ill)    = vo.spmcup    {ic}{iitt}(:,ill);    
%          vf.spmcdn    {ic}{iitt}(:,ill)    = vo.spmcdn    {ic}{iitt}(:,ill);    
%          vf.thetae    {ic}{iitt}(:,ill)    = vo.thetae    {ic}{iitt}(:,ill);  
%          vf.relhum    {ic}{iitt}(:,ill)    = vo.relhum    {ic}{iitt}(:,ill);  
%          vf.B         {ic}{iitt}(:,ill)    = vo.B         {ic}{iitt}(:,ill);
%          vf.enst      {ic}{iitt}(ill)      = vo.enst      {ic}{iitt}(ill);    
%          vf.frac_prec0{ic}{iitt}(ill)      = vo.frac_prec{ic}{iitt}(ill,1);
%          vf.frac_prec3{ic}{iitt}(ill)      = vo.frac_prec{ic}{iitt}(ill,2);
%          vf.frac_prec5{ic}{iitt}(ill)      = vo.frac_prec{ic}{iitt}(ill,3);
%else % glue new variables
%          prate = vo.prate     {ic}{iitt}(:,ill);
%          ctype = vo.ctype     {ic}{iitt}(:,ill); 
      %    sumprate = nansum(prate(mcsllx));
      %    vf.frac_precc{ic}{iitt}(ill) = nansum(prate(intersect(find(ctype==1),mcsllx)))/sumprate;
      %    vf.frac_precl{ic}{iitt}(ill) = nansum(prate(intersect(find(ctype==0),mcsllx)))/sumprate;
      %%    vf.frac_prec0{ic}{iitt}(ill)= nansum(prate(intersect(find(prate>0 & prate<3),mcsllx)))/numel(mcsllx);
      %%    vf.frac_prec3{ic}{iitt}(ill)= nansum(prate(intersect(find(prate>=3),mcsllx)))/numel(mcsllx);
%          vf.frac_prec0{ic}{iitt}(ill)= sum(prate(mcsllx)>=0.3)/numel(mcsllx);
%          vf.frac_prec3{ic}{iitt}(ill)= sum(prate(mcsllx)>=3)/numel(mcsllx);
%          vf.frac_prec5{ic}{iitt}(ill)= sum(prate(mcsllx)>=5)/numel(mcsllx);
%          vf.du03km_cam_abs{ic}{iitt}(ill) = vo.du03km_cam_abs{ic}{iitt}(ill);
%end
%  else 
%   vf.mcsillt4Cl{ic}{iitt}(ill,1:3) = vo.mcsillt4Cl{ic}{iitt}(ill,1:3);
%   vf.mcsllx4Cl{ic}{iitt}{ill} =      vo.mcsllx4Cl{ic}{iitt}{ill}  ;
%   vf.mpdt{ic}{iitt}(:,ill)     =     vo.mpdt{ic}{iitt}(:,ill)     ;
%   vf.mpdq{ic}{iitt}(:,ill)     =     vo.mpdq{ic}{iitt}(:,ill)     ;
%        vf.w{ic}{iitt}(:,ill)  = vo.w{ic}{iitt}(:,ill);
%        vf.theta     {ic}{iitt}(:,ill) = vo.theta     {ic}{iitt}(:,ill);
%        vf.mcsillt4Cl{ic}{iitt}(ill,1:3)  = vo.mcsillt4Cl{ic}{iitt}(ill,1:3);
%{
        vf.u_cam{ic}{iitt}(:,ill)  = vo.u_cam{ic}{iitt}(:,ill);
        vf.v_cam{ic}{iitt}(:,ill)  = vo.v_cam{ic}{iitt}(:,ill);
        vf.theta0{ic} = mean([vo.theta{ic}{:}],2); % mean profile over MCS lifecycle
        vf.thetaa{ic}{iitt}(:,ill) = vo.theta{ic}{iitt}(:,ill)-vf.theta0{ic}; % anomalous theta 
        vf.rho0{ic} = mean([vo.rho{ic}{:}],2); % mean profile over MCS lifecycle
        vf.thetaz0{ic} = mean([vo.thetaz{ic}{:}],2); 
%}
%        vf.pv     {ic}{iitt}(:,ill) = vo.pv{ic}{iitt}(:,ill)./vo.rho{ic}{iitt}(1:end-1,ill); % total pv
%        vf.pv0    {ic}{iitt}(:,ill) = parm.co(vo.mcsillt4Cl{ic}{iitt}(ill,2))*vf.thetaz0{ic}./vf.rho0{ic}(1:end-1); % pv MCS mean
%        vf.pva    {ic}{iitt}(:,ill) = vf.pv{ic}{iitt}(:,ill) - vf.pv0{ic}{iitt}(:,ill); % pv anomoly from MCS mean
  end 
