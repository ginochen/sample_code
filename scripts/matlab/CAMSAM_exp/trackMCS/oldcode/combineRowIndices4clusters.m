function [llcri, lonlatcri, lonlatsri] = combineRowIndices4clusters(llcri_it,lonlatcri_it,lonlatsri_it,nt)
  % Purpose: combine the timestep-separated mcs.llcentroid row indices into a single matrix to track a cluster
  % llcri_it{it}(icluster,1:2): row indices of mcs.llcentroid for each cluster (2-tsteps) 
  % llcri(icluster,it): combined matrix of llcentroid row indices for each cluster (timeseries)
  % lonlatcri(icluster,it): combined matrix of mcs.lonlat centroid row indices for each cluster (timeseries)
  % nt: # of total timesteps, e.g., [11:10:401]
  iic=1; % row counter for clusters
  for it = 1:nt-2 % only nt-1 llcri_it's, hence nt-1-1
    if it==1 % first loop
      ind = llcri_it{1}(:,1); % use all row indices for the 1st-tstep 
    else % after first loop
      ind = setdiff([llcri_it{it}(:,1);0],llcri_it{it-1}(:,2))'; % starting a new cluster if the rows aren't used previously
    end
    for iic = 1:numel(ind)
      i = ind(iic); % index for all it=1 clusters
     
      if (i & llcri_it{it}(i,2)~=0 ) % if two time centroids exists 
            llcri(iic,it:it+1) =     llcri_it{it}(i,1:2); % save the two centroid llcentroid-row to llcri matrix
        lonlatcri(iic,it:it+1) = lonlatcri_it{it}(i,1:2); % save the two centroid lonlat-row to lonlatcri matrix
        lonlatsri{iic,it}      = lonlatsri_it{i,it}{1};
        lonlatsri{iic,it+1}    = lonlatsri_it{i,it}{2};
        itt=1;
        while ( llcri_it{it+itt}(llcri(iic,it+itt),2)~=0) % if eq zero then go to the next i in ind
          ii = llcri(iic,it+itt); % llcentroid-row
          llcri(iic,it+itt+1)     =     llcri_it{it+itt}(ii,2);
          lonlatcri(iic,it+itt+1) = lonlatcri_it{it+itt}(ii,2);
          lonlatsri{iic,it+itt+1} = lonlatsri_it{ii,it+itt}{2};
          itt=itt+1; 
          if (it+itt==nt); break; end
        end
      end
    end
  end



!!!!!!!!!!!!!!!something's wrong with ii
and why is the 2nd column of llcri_it repeating!?
