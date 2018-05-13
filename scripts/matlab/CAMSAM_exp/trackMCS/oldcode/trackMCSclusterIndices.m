%%%% MCS / Cloud Cluster Centroid Tracking Algo %%%%%%%%%%%%%%
% Purpose: Obtain timeseries row indices (icr and isr) for the centroids of the
%          same clusters. Use the row indices to track a cluster's evolution
% Caution: Apply this code when varCAM0D-1D-2D variables are the same # of
%          lon-lat points as mcs.lonlat
%%%%%%%% PARAMETERS %%%%%%%%%%%%%
cdt = 10; % track cluster timestep gap
mint = 1; % minimum timesteps to qualify as mcs 3 steps = 3*5hrs = 15hrs 
maxd = 5; % max dist of two centroids to remain qualified as the same cluster (between two timesteps)
% according to Corfidi 2003, max propagation speed is around 60 knot ~ 30 m/s ~ 108 km/hr
t = [11:cdt:401]; %:10:401
nt= numel(t);
if (nt <= mint); error(['specify more than ' num2str(mint) ' steps!']); end

[llcri_it, lonlatcri_it, lonlatsri_it] = getRowIndices4Clusters(t,nt,maxd);
[llcri, lonlatcri, lonlatsri]       = combineRowIndices4clusters(llcri_it,lonlatcri_it,lonlatsri_it,nt);

%%%%%%%%% Get Var with Row Indices %%%%%%%%%%%%%%%%
clt = sum((llcri~=0),2); % # of times in a cluster 
cli = find(clt>=mint)'; % index of clusters that are greater than mint times 
ntClust = clt(cli); % the clt that is actually associated to the saved varCAM0D_mcs
if isempty(cli); error(['no clusters are greater than steps = ' num2str(mint)]); end
ii=1;
[varCAM0D_mcs, varCRM0D_mcs, varCAM1D_mcs, varCRM1D_mcs, varCRM2D_mcs, lonlatmcs] = getClusterVar(cli,ntClust,lonlatsri,lonlatcri,llcri,t);

save('var_PC1_MCSt','varCAM0D_mcs','varCAM1D_mcs',...
     'varCRM0D_mcs','varCRM1D_mcs','varCRM2D_mcs',...
     't', 'maxd', 'mint','llcri','llcri_it',...
     'lonlatsri','lonlatsri_it','lonlatcri','lonlatcri_it','cli','clt','ntClust','lonlatmcs')


