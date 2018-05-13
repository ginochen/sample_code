function mcs_corr(season,ilndocn)
% calc the correlation between MCS org parameters and duration or size of cluster

vname = 'dU03km_cam';
%vname = 'du03km_cam_abs';
%vname = 'du03km_abs';
%vname = 'freq_precc';
%vname = 'freq_precl';
%vname = 'prect'; % good corr
%vname = 'precc'; % good corr
%vname = 'Bave'; % no good corr
%vname = 'LIave'; % no corr
%vname = 'LImin'; % no corr
%vname = 'Bmax'; % no good corr
%ilndocn = 0;
%season = 'DJF';
%season = 'JJA';
it = 20; isize=5; itime=12;
%it = 10; isize=1; itime=6;
ilndocnstr = {'ocn','lnd'};
% obtain the clusters of nt4Cl that are non-end clusters (dropped in the code mcs_cluster_var)
eval(sprintf(['load /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/' season '/mcs_cluster_stats/mcs_stats_' ilndocnstr{ilndocn+1} ' nlndocnCl iclndocn ' vname ]))
% dU03km_cam %B LI du03km_abs 
ii=1;
maxnt = 20
for nt2=3:maxnt
  if nlndocnCl{nt2} <= 1; continue; end % cannot calc stats with one sample
  for iic = 1:nlndocnCl{nt2}
    icc(ii) = iclndocn{nt2}(iic); % cluster index after removing end-time cluster indices
    ii=ii+1;
  end
end

eval(sprintf(['load /Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly2/atm/hist/lat2525/' season '/mcs_clusters']));
%clear iCends
%for ic = 1:nCl;if any(ismember(t4Cl(ic,:),[1 numel(t)])); iCends(ic) = 1;else; iCends(ic)=0; end; end
%intersect(icc,find(iCends))
% calc corr of duration and MCS org parameters
%eval(sprintf(['corr(nt4Cl(icc)'',' vname '.mcscomposite.data(it,:)'')']))
%tmp = nt4Cl(island==0 & ~iCends')+2 <=maxnt;
corr(nt4Cl(icc)',mcsnllave(icc)')

figure
eval(sprintf(['vtmp = (' vname '.mcscomposite.data);']))
%corr((nt4Cl(icc)'+2)*3,nanmean(vtmp,1)') % corr between duration and averaged strenght of variable
%corr((mcsnllave(icc)'+2)*3,nanmean(vtmp,1)') % corr between duration and averaged strenght of variable
%eval(sprintf(['plot(nanmean(' vname '.mcscomposite.data,2))'])); % average LIave over the ocean or land

figure
tmp = (nt4Cl(icc)+2)*3; 
i = tmp>itime ;disp([num2str(sum(i)) ' of time > ' num2str(itime) ' hours clusters'])
eval(sprintf(['[cv pv] = corr(tmp(i)'',nanmean(' vname '.mcscomposite.data(:,i),1)'',''type'',''pearson'',''rows'',''complete'')']))
%eval(sprintf(['corr(tmp(i)'',' vname '.mcscomposite.data(it,i)'',''type'',''pearson'')']))
%eval(sprintf(['corr(tmp(i)'',max(abs(' vname '.mcscomposite.data(:,i)))'',''type'',''pearson'')']))
eval(sprintf(['scatter(tmp(i)'',' vname '.mcscomposite.data(it,i)'')']))

% check the percentage of long-live cell has shear greater than a threshold
% long lived cell cluster index = i
eval(sprintf(['for ii=1:20; nsh(ii) = sum(' vname '.mcscomposite.data(ii,i)>5); end'])) % set threshold manually here
disp(['percentage of ' vname ' greater than a threshold ' num2str(100*mean(nsh/sum(i))) '%']); % percentage is slightly greater over land

figure
% calc corr of size and MCS org parameters
%eval(sprintf(['corr(mcsnllave(icc)'',' vname '.mcscomposite.data(it,:)'')']))
tmp = mcsnllave(icc);
i = i & tmp<16; %tmp>isize & tmp<14; 
disp([num2str(sum(i)) ' of size > ' num2str(isize) ' deg^2 clusters'])
eval(sprintf(['corr(tmp(i)'',nanmean(' vname '.mcscomposite.data(:,i),1)'',''type'',''pearson'')']))
%eval(sprintf(['corr(tmp(i)'',' vname '.mcscomposite.data(it,i)'',''type'',''pearson'')']))
%eval(sprintf(['corr(tmp(i)'',max(abs(' vname '.mcscomposite.data(:,i)))'',''type'',''pearson'')']))
eval(sprintf(['scatter(tmp(i)'',' vname '.mcscomposite.data(it,i)'')']))
title('size')

% for cluster size larger than 6 deg^2 corr = 0.2575
% for cluster time longer than 10 hours corr = 0.19
