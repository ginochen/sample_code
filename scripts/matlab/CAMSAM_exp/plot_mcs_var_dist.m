function var = plot_mcs_var_dist(vname,season,dosave,donewplot)
% output var to calc corr
  docbar = 1;
  if ~exist('donewplot')
    donewplot=1; % hold on old plot (0) or open a new figure (1)
  end
  if ~exist('dosave')
    dosave = 0
  end
  figure('units','normalized','outerposition',[0 0 1 1])
  lndocnstr = {'ocn','lnd','lndocn','nocond'};
  casei='F_2000_SPCAM_m2005_3hrly2'
%  casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch'
%  diri = ['/Users/g/archive/matlab/figure/var/mcsVars/mcs_3hrly/f09f09/'];
%  diro = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/' season];
  dlat = 'lat2525';
  diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season]
  diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/'] % output dir for figures
  dircbar = ['/Users/g/archive/matlab/figure/mcs_cluster/f09f09/'];
  load([diri '/mcs_clusters.mat'],'nt4Cl','t','t4Cl','mcsnllave','mcsnll','nCl','lon','lat','ilatlim','latlim');
  load([diri '/mcs_cluster_parm.mat']);
  load([diri '/mcs_cluster_var/mcs_cluster_' vname '.mat'])
  eval(['vtmp =' vname ';']);
%  load([diri '/mcs_cluster_var/mcs_cluster_ctype.mat'])
%  load([diri '/mcs_cluster_var/mcs_cluster_frac_precc.mat'])
%  load([diri  '/mcs_cluster_var/mcs_cluster_frac_precl.mat'])
%  load([diri  '/mcs_cluster_var/mcs_cluster_precc.mat'])
%  load([diri  '/mcs_cluster_var/mcs_cluster_precl.mat'])
  load([diri '/mcs_cluster_var/mcs_cluster_mcsillt4Cl.mat'])
  vnames = {'du03km_abs','du03km_cam_abs','dU03km_cam','precc','precl','prect','LIave','enst','fftke','Bave'}; 
%  cbarunits={'[m s^{-1}]','[mm hr^{-1}]','[mm hr^{-1}]'};
  cbarunits={'[m s^{-1}]','[m s^{-1}]','[m s^{-1}]','[mm]','[mm]','[mm]','[K]','[m s^{-1}]','[J kg^{-1}]','[m s^{-1}]'};
  ncon = 200;
  icmap = [25,25,25,12,12,12,12,25,25,25];
  doflip = [0,0,0,0,0,0,1,0,0,0];
  vdim = containers.Map([vnames],[zeros(1,numel(vnames))]);
  ivv  = containers.Map([vnames],[1:numel(vnames)]);
  caxislim={[0 10],[0 10],[0 10],[0 15],[0 15],[0 7],[-5 0],[0 3e-2],[0 2],[0 40]};
  if ismember(vname,{'fftke'})
    ixcoef=2;
    km2m=1000;
    xwavelength = parm.nx*parm.dx./(1:parm.nx/2)/km2m; % a_0 + b_0 + \sum ( a_n * cos(2pi/(L/n)*x) + b_n * sin(2pi/(L/n)*x) ), where L=nx*dx, n=1:nx/2
    zwavelength = parm.nzi*parm.dzi./(1:parm.nzi/2)/km2m; % the entire height is 28km, first wave covers entire 28km with one wavelength, so zwavelength=nzi*dzi
    freq = 1./xwavelength; % [km^-1]
    if ixcoef>1
      disp(['wavelength = ' num2str(xwavelength(ixcoef-1))])
      disp(['frequency = ' num2str(freq(ixcoef-1)) '[km^{-1}]'])
    end
    x = [1:parm.nx/2] - 1; % -1 to start at zero
  end
  %varmap = zeros(288,192);
  %countmap = zeros(288,192);
  varmap = zeros(288,192);
  countmap = zeros(288,192);
  nt = numel(t);
  ii=1;
  for ic=1:nCl
    if any(ismember(t4Cl(ic,:),[1 nt])); continue; end % skip the cluster that starts/ends at the first/last time step since cannot index the -1/+1 time step 
    for iitt=1:numel(mcsillt4Cl{ic})%nt4Cl(ic)+2
      for ill=1:size(mcsillt4Cl{ic}{iitt},1)
        l1 = mcsillt4Cl{ic}{iitt}(ill,1);
        l2 = mcsillt4Cl{ic}{iitt}(ill,2);
        if ~ismember(vname,{'fftke'})
          varmap(l1,l2) = vtmp{ic}{iitt}(ill)+varmap(l1,l2);
          if numel(mcsillt4Cl{ic})>3 % greater than 3 timesteps = 9 hours
            var(ii) = vtmp{ic}{iitt}(ill); ii=ii+1; % this is outputted to calc corr
          end
        else
          varmap(l1,l2) = vtmp{ic}{iitt}(ixcoef,1,ill)+varmap(l1,l2); % largest wave at zonal wave# 1 and vertical wave# 0
        end
        countmap(l1,l2) = countmap(l1,l2)+1;
      end
    end
  end
  disp(['total mcs = ' num2str(sum(countmap(:)))])
%  if strcmp(stats,'mean')
  if ~ismember(vname,{'precc','precl'})
    idx = [countmap>0]; % total # of non-zero idx doesn't need to be nCl, can be less or more
    varmap(idx) = varmap(idx)./countmap(idx);
  else % precc precl
    varmap = varmap/2; % divide by two for mm/30min, and added for seasonal 3 mm/30min * 30min + 1.5 mm/30min * 30min + ... + 1.2 mm/30min * 30min = X mm
  end
%  elseif strcmp(stats,'median')
  if donewplot
    figure('units','normalized','outerposition',[0 0 1 1])
  else
    hold on
  end
%  if ismember(vname,{'enst'})
%    contourf(lon,lat(ilatlim(1):ilatlim(2)),varmap(:,ilatlim(1):ilatlim(2))',ncon);%,[0.001:0.005:0.1],'color',[0.75 0.75 0.75]);
%  else 
  varmap(varmap==0)=NaN; 
    contourf(lon,lat(ilatlim(1):ilatlim(2)),varmap(:,ilatlim(1):ilatlim(2))',ncon,'linestyle','none');
%  end
  colormap(cmap(icmap(ivv(vname)),doflip(ivv(vname)))); hold on
  ylim(latlim);
  title([season])
  coast_centered(0); % center the coastline at 0 degrees
%  caxislim=[0,150];
  caxis(caxislim{ivv(vname)});
  set(gca,'Fontsize',30); 
  daspect([1 .5 1])
  if ismember(vname,{'fftke'})
    vnameo = ['fftke_' num2str(ixcoef)];  
  else
    vnameo = vname;
  end
  if dosave
%    saveas(gcf,[diro '/mcs_' vnameo '_dist.fig'])
    disp(['saveas(gcf,' diro '/mcs_' vnameo '_dist.png)'])
    saveas(gcf,[diro '/mcs_' vnameo '_dist.png'])
    crop([diro '/mcs_' vnameo '_dist.png'])
  end
% COLORBAR
  if docbar
    figure('units','normalized','outerposition',[0 0 1 1])
    h = colorbar('southoutside'); 
    axis off
    h.Label.String = cbarunits{ivv(vname)}; ;
    h.Position = [0.05 0.1170 0.7750 0.0181];
    h.Label.Position = [caxislim{ivv(vname)}(2)+0.1*abs(diff(caxislim{ivv(vname)})) 0 0];
    set(gca,'Fontsize',44);
%    set(get(h,'title'),'string','[Cells]')
%    cmapp = cmap(icmap(ivv(vname)));
%    icmapp = round(linspace(1,size(cmapp,1),ncon));
%    colormap(cmapp(icmapp,:));
    colormap(cmap(icmap(ivv(vname)),doflip(ivv(vname))));
    caxis(caxislim{ivv(vname)});
%    saveas(gca,[dircbar '/mcs_' vnameo '_dist_colorbar.fig'])
    saveas(gca,[dircbar '/mcs_' vnameo '_dist_colorbar.png'])
          crop([dircbar '/mcs_' vnameo '_dist_colorbar.png'])
  end
%
  disp('dosave=0; vname=''precc''; season={''JJA'',''DJF''}; for is = 1:numel(season); plot_mcs_var_dist(vname,season{is},dosave); end')
  disp('dosave=0; vname=''du03km_abs''; season={''JJA'',''DJF''}; for is = 1:numel(season);plot_mcs_var_dist(vname,season{is},dosave); end')
  disp('dosave=0; vname=''B''; season={''JJA'',''DJF''}; for is = 1:numel(season);plot_mcs_var_dist(vname,season{is},dosave); end')
