% Purpose: Calculate variable on PC1 lon-lat indices  
%
run ~/scripts/matlab/startup.m
load(['/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/var_PC1_' num2str(it) '.mat'],'var','div','vort','FFTke','parm');
nx=32;
ivec = 1; % 
idx.u=1;
for ib=1:4
   iCAMll{ib} = find(parm.iillt_zm_PC{ib}(3,:)==it);  % find the it timestep large-PC column indices
   iCAMllNH{ib} = find(lat(parm.iillt_zm_PC{ib}(2,iCAMll{ib}))>0); % find the NH column indices in iCAMll{ib} to ammend vorticity sign issue in SH
   iCAMllSH{ib} = find(lat(parm.iillt_zm_PC{ib}(2,iCAMll{ib}))<0); % find the SH column indices 
   ii=1; % counter for saved PC1 indices
   for ill = iCAMll{ib}
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Calc lag-lead U(x,z) field EOF1 for all lon-lat indices and take a %
      % spatial average to get the composite                               %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      u_tmp = reshape(var{ib}(1:nx,1:parm.nzi-5,1:parm.nlag*2+1,ii,idx.u),nx*(parm.nzi-5),parm.nlag*2+1);
      [U S V] = svd(u_tmp');
      % EOF of u field over the lag-lead time for one lat-lon point
      u_tlag_eof{ib}(:,:,ii) = reshape(mean(U(:,ivec)*S(ivec,ivec))*V(:,ivec),nx,parm.nzi-5); % mean U at eof_$ivec over the nlag-nlead time
      ii=ii+1;
   end
   %%%%%%%%%%%%%%%%%%%%%%
   % k-means clustering %
   %%%%%%%%%%%%%%%%%%%%%%
   u_tmp = reshape(u_tlag_eof{ib},nx*(parm.nzi-5),size(u_tlag_eof{ib},3)); % matrix_in needs to be nsamples-by-pdim
%   nk=10; % separate into nk clusters
%   [idx C] = kmeans(u_tmp',nk); %idx: predicted cluster indices; C: k-by-p centroid location
%   for ik=1:nk;
%      ik
%      contourf(reshape(C(ik,:),nx,parm.nzi-5)');
%      pause
%   end
   %%%%%%%%%%%%%%%%%%%%%%
   %  90% variance eof  %
   %%%%%%%%%%%%%%%%%%%%%%
   [svec sval pc] = svdVar(u_tmp);
   pc_ave = mean(pc,1); % average over all lon-lat indices (originally ave over time)
   for ieof = 1:size(pc,2)
      ieof
      u_tlag_eof_composite{ib}(1:nx,1:parm.nzi-5,ieof) = reshape(pc_ave(ieof)*svec(:,ieof),nx,parm.nzi-5);
      contourf(u_tlag_eof_composite{ib}(:,:,ieof)'); pause
   end
end
save('Basinwise-omposite','u_tlag_eof_composite','u_tlag_eof','idx_NH','idx_SH','iCAMll');
['done calc u_tlag_eof for it=' num2str(it)]
