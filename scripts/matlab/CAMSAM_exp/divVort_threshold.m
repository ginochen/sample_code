% Purpose: find the column indices of iillt_zm_PC{ib} that has threshold_bigPC >= 20 
%
% Ex. find the variables that are 20+ by
% vort{ib}(:,:,:,iCAMll_thres{ib}(:,ithr))
% 
% Ex. for threshold 20 at 1-th index, we can find the lon-lat-time
% index set by:
% iillt_zm_PC{ib}(:,nonzeros(iCAMll{ib}(iCAMll_thres{ib}(:,1))))
%
%
load '/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat'
%load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/PC1Threshold10/var_PC1_11.mat','parm')
thres_bigPC = [20:10:100]; % thres_bigPC*0.5 will give the maximum K/day heating rate in the vertical
iEOF=1; 
bsign=[1 -1 1 1];
tvec = [21:10:491];
for it = tvec
   tic
   load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/var_PC1_' num2str(it) '.mat'],'parm')
   iillt_zm_PC = parm.iillt_zm_PC;
   for ib=1:4
      % find the columns in iillt_zm_PC that are first at 'it', second threshold 20+
      iCAMll{ib} = find(iillt_zm_PC{ib}(3,:)==it);  % find columns of iillt_zm_PC that are at it-timestep
      tt = 1; % counter for number of indices in iCAMll exceeding all 20+ thresholds 
      for ithr = 1:length(thres_bigPC)
         ii=1; % counter for number of indices in CAMll exceeding ithr-threshold 
         for i=1:length(iCAMll{ib}) 
            icol = find(iillt_zm{ib}(1,:)==iillt_zm_PC{ib}(1,iCAMll{ib}(i)) & ...
                        iillt_zm{ib}(2,:)==iillt_zm_PC{ib}(2,iCAMll{ib}(i)) & ...
                        iillt_zm{ib}(3,:)==iillt_zm_PC{ib}(3,iCAMll{ib}(i))); 
            if ( bsign(ib)*pc_spz{ib}(icol,iEOF) > thres_bigPC(ithr) )
               ni(tt) = ii; nj(tt) = ithr; nv(tt) = i;
               ii=ii+1; tt=tt+1;
            end
         end
      end
      iCAMll_thres{ib} = sparse(ni,nj,nv);  % the index of iCAMll{ib} that are 20+
      clear ni nj nv
   end
   parm.iCAMll = iCAMll;
   parm.iCAMll_thres = iCAMll_thres;
   parm.thres_bigPC = thres_bigPC;
   parm.bsign = bsign;
   save(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/var_PC1_' num2str(it) '.mat'],'parm','-append')
   toc
end
