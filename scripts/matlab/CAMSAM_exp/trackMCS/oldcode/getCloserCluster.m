
function ic2 = getCloserCluster(ic2,isr1,mcs1,mcs2); 
  % Purpose: get the closer cluster to mcs1 if there are two centroids in mcs2 that are the same min dist
  % Method: use 'pdist2' to eval the mean Euclidean dist of all surrounding points of mcs1 and mcs2
  % ic2: llcentroid row index of 2
  for iii=1:numel(ic2)
    [icr2 isr2] = getRowIndices4Cluster(ic2(iii),mcs2);
    pdist_pts = pdist2(mcs1.lonlat(isr1,:),mcs2.lonlat(isr2,:));
    meanDist(iii) =mean(mean(pdist_pts)); % mean distance between all points of the two time MCS 
    overlapPt(iii) = numel(intersect(isr1,isr2)); % overlapped points
  end 
  if any(overlapPt) | any(pdist_pts<1.6) % overlapping points exists or any point within 1 grid point distance
    ic2 = ic2(find(meanDist == min(meanDist) & overlapPt==max(overlapPt)));% eval which mean distance is closer & has the most overlapped points
    if numel(ic2)>1 % if still multiple clusters satisfies, then use the first one arbitrarily for now
      ic2 = ic2(1);
    end
  else
    ic2=0;
  end
