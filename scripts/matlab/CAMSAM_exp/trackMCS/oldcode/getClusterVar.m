

function [varCAM0D_mcs, varCRM0D_mcs, varCAM1D_mcs, varCRM1D_mcs, varCRM2D_mcs, lonlatmcs] = getClusterVar(cli,ntClust,lonlatsri,lonlatcri,llcri,t)
  ii=1;
  for ic = 1:numel(cli)
    iit = find(llcri(cli(ic),:)~=0); % llcri(icluster,it) saves the mcs.llcentroid row indices as non-zero values, so finding the non-zero columns of llcri is finding the timesteps assoc to a cluster 
    for it = 1:ntClust(ic)
      eval(sprintf('load(''var_PC1_%d.mat'',''varCAM0D'',''varCRM0D'',''varCAM1D'',''varCRM1D'',''varCRM2D'',''mcs'');',t(iit(it))));
      illr = lonlatsri{cli(ic),iit(it)}; % lonlat row indices in a cluster at a timestep
      lonlatmcs.centroid{ic}(:,it) = mcs.lonlat(lonlatcri(cli(ic),iit(it)),:); % save the centroid lon lat time-series for each cluster separately
      lonlatmcs.npts(ic,it) = numel(illr); % # of lon-lat points at a given tstep for a cluster
      for is = 1:lonlatmcs.npts(ic,it) % surroundings' row indices counter
        ir = illr(is); % lonlat row for surroundings' 
        % Analyze at every cluster/timestep 
        varCAM0D_mcs{ic}{it}(:,is)          = varCAM0D(ir,:);
        varCRM0D_mcs{ic}{it}(:,is)          = varCRM0D(ir,:);
        varCAM1D_mcs{ic}{it}(:,:,is)        = varCAM1D(:,ir,:);
        varCRM1D_mcs{ic}{it}(:,:,is)        = varCRM1D(:,ir,:);
        varCRM2D_mcs{ic}{it}(:,:,:,is)      = varCRM2D(:,:,ir,:);
        lonlatmcs.surrounding{ic}{it}(:,is) = mcs.lonlat(ir,:); % save the surrounding lon lat time-series for each cluster separately
        ii=ii+1;
      end
    end
  end
end
