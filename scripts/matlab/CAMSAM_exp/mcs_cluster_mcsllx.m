season = 'JJA'; 
  spcase      = 'F_2000_SPCAM_m2005_3hrly1';
  spArchive   = ['/projects/rsmas/kirtman/gchen/cesm_spcam/archive/'];
  spHist.dir  = [spArchive spcase '/atm/hist/'];
load([spHist.dir season '/mcs_clusters.mat'])
diro = [spHist.dir season '/mcs_cluster_var/'];
nt = numel(t);
tic
  for it=2:numel(t)-1
    for ic = 1:nCl
      if any(ismember(t4Cl(ic,:),[1 nt])); continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
      if (it == mcsillt4Cl{ic}{1}(1,3)) % cluster's starting time index
        for ill=1:mcsnll(ic,1) % go through the ilon-ilat at iit time for cluster ic
%  vo.mcsllx4Cl{ic}{iitt}{ill} = mcsllx4Cl{ic}{end}{ill};
          vf.mcsllx4Cl{ic}{1}{ill} = mcsllx4Cl{ic}{1}{ill};
          vf.mcsllx4Cl{ic}{2}{ill} = mcsllx4Cl{ic}{1}{ill};
          if (it == mcsillt4Cl{ic}{end}(1,3)) % for short lived clusters with only one time index 
            vf.mcsllx4Cl{ic}{3}{ill} = mcsllx4Cl{ic}{end}{ill};
          end
        end
      elseif (it < mcsillt4Cl{ic}{end}(1,3) & it > mcsillt4Cl{ic}{1}(1,3))% it's in between start and end of a cluster
        itt = find(t4Cl(ic,:)==it);
        for ill=1:mcsnll(ic,itt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
          vf.mcsllx4Cl{ic}{itt+1}{ill} = mcsllx4Cl{ic}{itt}{ill};
        end
      elseif (it == mcsillt4Cl{ic}{end}(1,3)) % cluster's ending time index
        ntt = nt4Cl(ic);
        for ill=1:mcsnll(ic,ntt) % sum(mcsnll(ic,:)~=0) is the last nonzero index of mcsnll(ic,:) 
          vf.mcsllx4Cl{ic}{ntt+1}{ill} = mcsllx4Cl{ic}{end}{ill};
          vf.mcsllx4Cl{ic}{ntt+2}{ill} = mcsllx4Cl{ic}{end}{ill};
        end
      end % if loop of calc variables 
    end
  end
toc
mcsllx4Cl = vf.mcsllx4Cl; 
save([diro '/mcs_cluster_mcsllx4Cl.mat'],'mcsllx4Cl')
