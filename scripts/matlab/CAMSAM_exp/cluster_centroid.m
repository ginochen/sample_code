function cluster_centroid(fout,season)
%%%%% Find cluster centroid time series %%%%%%%%%%%%% 
%  fo = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/F_2000_SPCAM_m2005_3hrly1/atm/hist/' season '/mcs_clusters.mat'];
%  res = 'f09f09';
  %fo = ['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/' res '/' season '/mcs_clusters.mat'];
  load(fout);
  display('start finding cluster centroids')
  % index mcsllt clusters by mcsllt(nonzeros(mcslltClRowInd(ic,:)),1:3)
  disp(nCl)
  for ic = 1:nCl
    for iit = 1:nt4Cl(ic)
      if size(mcsillt4Cl{ic}{iit},1)>1
        [~, mcsilltcentroids{ic}(iit,1:2)] = kmedoids(mcsillt4Cl{ic}{iit}(:,1:2),1,'algorithm','small','distance','euclidean'); % save the centroid lon-lat-time
        [~,illtcentroid{ic}(iit)] = ismember(mcsilltcentroids{ic}(iit,1:2),mcsillt4Cl{ic}{iit}(:,1:2),'rows'); % the row associated to centroid MCS in a ic-cluster
      else
        mcsilltcentroids{ic}(iit,1:2) = mcsillt4Cl{ic}{iit}(:,1:2);
        illtcentroid{ic}(iit) = 1;
      end
      mcsilltcentroids{ic}(iit,3) = t4Cl(ic,iit);
    end
  end
  display('end finding cluster centroids')
  % save file, and ask if -append variable to existing file
  save(fout,'mcsilltcentroids','illtcentroid','-append')
