load('var_PC1_11.mat','parm') % load a general parm to get the variable index 'idv' and 'parm.zint'
load var_PC1_MCSt.mat
idv = getVarIndex({'varCRM2D','varCAM1D'},parm);
%nmode = size(svec{ib},2); % all modes
nbin = 200;
nmode = 3; % # of modes to plot
zint = parm.zint/1000; % height coordinate 
nr = 1; % # of rows for subplot
nc = nmode; % # of columns for subplot
nt = mint+1; % # of timesteps in a cluster
tvec = find(ntClust>=nt); 
%%%%% REUSE varCAM*D_mcs FOR PLOTTING %%%%%%%%%%%%%
var1D_mcsall = selectLargeTimeCluster(varCAM1D_mcs, ntClust, nt); % Lump all cluster/timestep into one big matrix
figure; 
plot_CAFD(squeeze((var1D_mcsall(:,idv.varCAM1D.dtcond,:))), [nr nc], 1,  [-10 -5 -5], [10 5 5], nbin, nmode, zint); % CAM subplots
figure;
plot_CAFD(squeeze((var1D_mcsall(:,idv.varCAM1D.spdt,:))), [nr nc], 1, [-10 -5 -5], [10 5 5], nbin, nmode, zint); % SPCAM subplots
%var1D_mcsall_big = selectLargeSizeCluster(varCAM1D_mcs, ntClust, 2, lonlatmcs); % Lump all cluster/timestep into one big matrix
%figure; 
%plot_CAFD(squeeze((var1D_mcsall_big(:,idv.varCAM1D.dtcond,:))), [nr nc], 1,  [-10 -5 -5], [10 5 5], nbin, nmode, zint); % CAM subplots
%figure; 
%plot_CAFD(squeeze((var1D_mcsall_big(:,idv.varCAM1D.spdt,:))), [nr nc], 1,  [-10 -5 -5], [10 5 5], nbin, nmode, zint); % CAM subplots




%%%%%%%%%%%% FUNCTIONS BEGIN %%%%%%%%%%%%%%%%%%%%%%%%%%

function [var_all] = selectLargeSizeCluster(var, ntClust, npt, lonlatmcs)
  % Purpose: Lump all cluster, timesteps, surrounding points (according to # of gridpoint threshold) into one big matrix for svd and scatter plot use
  numdim = ndims(var{1}{1});
  ii=1; % counter for the lumped index
  maxClustPts = max(lonlatmcs.npts,[],2); % maximum points at a time for each cluster
  for ic = 1:numel(ntClust)
    for it = 1:ntClust(ic)
      if (maxClustPts(ic) >= npt) % filter the maximum points covering over the size threshold
        for is = 1:size(var{ic}{it},numdim) % total number of surrounding points at a timestep for a cluster
          if (numdim==2) % varCAM0D_mcs, varCRM0D_mcs
            var_all(:,ii) = var{ic}{it}(:,is);
          elseif (numdim==3) % varCAM1D_mcs, varCRM1D_mcs
            var_all(:,:,ii) = var{ic}{it}(:,:,is);
          elseif (numdim==4) % varCRM2D_mcs
            var_all(:,:,:,ii) = var{ic}{it}(:,:,:,is);
          else
            error('dimension of var is wrong')
          end
          ii=ii+1;
        end
      end
    end
  end
end

function [var_all] = selectLargeTimeCluster(var, ntClust, nt)
  % Purpose: Lump all cluster, timesteps, surrounding points (according to life-time threshold) into one big matrix for svd and scatter plot use
  % varin:
  %   ntClust: total # of timesteps in a cluster
  %   nt:      minimum threhold timesteps for a cluster to be selected
  numdim = ndims(var{1}{1});
  ii=1; % counter for the lumped index
  for ic = 1:numel(ntClust)
    if ( ntClust(ic) >= nt ) % if cluster total life-timesteps greater than nt
      for it = 1:ntClust(ic)
        for is = 1:size(var{ic}{it},numdim) % total number of surrounding points at a timestep for a cluster
          if (numdim==2) % varCAM0D_mcs, varCRM0D_mcs
            var_all(:,ii) = var{ic}{it}(:,is);
          elseif (numdim==3) % varCAM1D_mcs, varCRM1D_mcs
            var_all(:,:,ii) = var{ic}{it}(:,:,is);
          elseif (numdim==4) % varCRM2D_mcs
            var_all(:,:,:,ii) = var{ic}{it}(:,:,:,is);
          else
            error('dimension of var is wrong')
          end
          ii=ii+1;
        end
      end
    end
  end
end

function plot_CAFD(variable, nsubp, ii, xmin, xmax, nbin, nmode, z)
  % ii:  subplot start row
  % xmin: min binedge for histogram
  % xmax: max binedge for histogram
  % nbin: # of bins to do histogram at each height level
  % z: height coordinate to plot projected EOFs
  % nsubp: [nr nc] # of rows and columns for subplot 
  [svec var pc pcmean] = svdVar(variable,1,1,0,1); 
  nz = size(svec,1);
  zvec = [1:nz]; % use the surface (nz) to dz+1-th (nz-dz) pressure level, notice pressure increases with index
  ntsamp = size(pc,1);
  fracvar = var/sum(var);
  for im=1:nmode
    xkb = pcmean(im)*svec(zvec,im);
    for iz = 1:nz
      xk_orig(iz,:) = xkb(iz) + pc(:,im)*svec(zvec(iz),im);  
    end
    subplot(nsubp(1),nsubp(2),ii)
    [vprob, vbincenters] = binVarWithHeight(xk_orig,xmax(im),xmin(im),nbin,nz);
    contourf(vbincenters,z,vprob,1000,'linestyle','none')
    hold on; plot(xkb,z,'r-')
    text(vbincenters(end-2),2,num2str(fracvar(im)))
    colormap(cmap(5))
    display(im)
    %display(['fractional variance = ' num2str(fracvar(im))])
    ii=ii+1;
  end
end
