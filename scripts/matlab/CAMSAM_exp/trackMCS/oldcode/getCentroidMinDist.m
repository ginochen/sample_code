function [dist mind] = getCentroidMinDist(mcs1,mcs2,ic1)
  % Purpose: get the minimum dist between mcs1 (ic1-th centroid) and mcs2 (all centroids)
  dist = sqrt((mcs1.llcentroids(ic1,1)-mcs2.llcentroids(:,1)).^2+...
         (mcs1.llcentroids(ic1,2)-mcs2.llcentroids(:,2)).^2);
  mind = min(dist);
