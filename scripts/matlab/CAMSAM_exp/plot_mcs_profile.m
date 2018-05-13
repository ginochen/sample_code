function [svec pc] = plot_mcs_profile(vname,stats,ilndocn,ilatzone,season,dlat,ifhold,ilife,dosave,ii)
% vname='spdt';ilatzone=3;dlat='lat9090';season='JJA';dosave=0;ii=1;for ilnd=[0 1]; for ilife=[1 10 20];eval(sprintf('[svec.%s{ii} pc.%s{ii}]=plot_mcs_profile(vname,    ''eof'',ilnd,ilatzone,season,dlat,0,ilife,dosave,ii);',vname,vname)); ii=ii+1;end;end;
% [R{ii} P{ii}] = corr(pc_spdq{ii}(:,1:nmode),pc_camdq{ii}(:,1:nmode));ii=ii+1;
% To do ddt or ddq, use vname={'spdt','camdt'}; vname={'spdq','camdq'}; otherwise set vname to single string, e.g., {'spdt'}
% [svec pc] = plot_mcs_profile(vname,stats,ilndocn,ilatzone,season,dlat,ifhold,ilife,signflip)
%function plot_mcs_profile(vname,stats,ilndocn,season,ifhold,ilife,signflip)
%function plot_mcs_profile(vname,season,stats)
%
% stats: 'eof', 'fcbt', 'fft'
% signflip: n-dim ones vector for n-mode, e.g., [1 -1 1] for 3-modes
% ilife: [1 20];
  if dosave
    dosevsigntest=0;
  else
    dosevsigntest=1;
  end
  dostatdata=0; %  use data (0)without (1)with spatial averaging and lifetime interp
  imode = 2;
%{
  if ismember(stats,{'eof'})
    ipm=17; % 17 for 15km, 20 for 18km
  elseif ismember(stats,{'fcbt'})
    ipm=20;
  end
%}
  lw = 5 ; % linewidth
  vnames = {'spdt','camdt','spdq','camdq'};
  scale = containers.Map(vnames,86400*[1,1,-2260,-2260]);
% sign flip manually according to the sign of pcmean*svec
if strcmp(vname,'ddt')
  vname1={'spdt','camdt'};
elseif strcmp(vname,'ddq')
  vname1={'spdq','camdq'};
else
  vname1={vname};
end
   
if strcmp(season,'DJF')  
sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
                                           {[-1,-1,1;-1,1,-1;1,-1,-1;-1,1,-1],...
                                            [1,-1,1;1,-1,-1;1,-1,-1;1,1,1],...
                                            [1,1,1;-1,-1,1;1,1,1;1,1,1],...
                                            [1,-1,-1;1,-1,-1;1,-1,1;1,-1,-1],...
                                            [-1,1,1;1,-1,1;-1,1,1;1,1,1],...
                                            [1,-1,1;1,-1,1;1,-1,-1;1,-1,1]});
elseif strcmp(season,'JJA')
%{
sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
                                           {[1,1,1; -1,1,1; -1,-1,-1;  -1,-1,1; 1,1,-1; -1,-1,-1],...
                                            [-1,-1,-1; 1,-1,-1; 1,-1,-1; 1,1,1; -1,1,1; -1,-1,-1],...
                                            [-1,-1,1; 1,-1,-1; 1,-1,-1; -1,1,-1; -1,1,1; -1,1,-1],...
                                            [1,1,1; 1,1,-1; 1,1,1; 1,1,1; -1,1,-1; -1,-1,1],...
                                            [-1,1,1; 1,1,1; 1,-1,1; -1,1,1; 1,1,1; 1,1,1],...
                                            [1,-1,1; 1,1,1; 1,-1,1; 1,-1,-1; 1,1,1; 1,-1,1]});
sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
                                           {[1,1; -1,1; -1,-1;  -1,-1; 1,1; -1,-1],...
                                            [-1,1; 1,1; 1,1; 1,-1; -1,-1; -1,1],...
                                            [-1,1; 1,1; 1,1; -1,-1; -1,-1; -1,-1],...
                                            [1,-1; 1,-1; 1,-1; 1,-1; -1,-1; 1,-1],...
                                            [-1,1; 1,1; 1,-1; -1,1; 1,1; 1,1],...
                                            [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
% all
%}
if ilatzone==3 % global
  if dostatdata
    sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
       {[-1,1; -1,1; -1,1; 1,1; 1,1; -1,1],...
        [1,1; 1,-1; 1,-1; -1,1; -1,1; -1,1],...
        [1,1; -1,1; -1,1; 1,1; -1,1; 1,-1],...
        [-1,1; 1,1; 1,1; 1,-1; -1,-1; 1,-1],...
        [-1,-1; -1,-1; -1,-1; -1,1; -1,-1; 1,1],...
        [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
  else
    sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
       {[-1,-1; 1,1; -1,-1; 1,1; 1,-1; 1,-1],...
        [1,1; 1,1; 1,-1; -1,-1; 1,-1; -1,-1],...
        [1,-1; 1,-1; -1,1; -1,1; -1,1; 1,1],...
        [-1,-1; 1,-1; 1,-1; -1,-1; -1,1; -1,-1],...
        [-1,-1; -1,-1; -1,-1; -1,1; -1,-1; 1,1],...
        [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
  end
% tropics
elseif ilatzone==0 
  if dostatdata
    sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
       {[1,1; -1,-1; 1,1; -1,1; -1,-1; -1,-1],...
        [1,-1; -1,1; 1,-1; -1,-1; 1,1; -1,1],...
        [1,-1; -1,-1; -1,-1; -1,1; -1,-1; -1,-1],...
        [-1,1; 1,1; 1,1; 1,-1; -1,-1; 1,-1],...
        [-1,-1; -1,-1; -1,-1; -1,1; -1,-1; 1,1],...
        [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
  else
    sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
       {[-1,-1; 1,1; -1,1; 1,1; 1,-1; 1,-1],...
        [1,-1; 1,-1; 1,-1; -1,-1; 1,1; -1,-1],...
        [1,-1; 1,1; -1,-1; 1,-1; -1,-1; 1,1],...
        [-1,1; 1,1; 1,1; 1,-1; -1,1; 1,1],...
        [-1,-1; -1,-1; -1,-1; -1,1; -1,-1; 1,1],...
        [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
  end
elseif ilatzone==1
    sign_vec = containers.Map({'spdt','spdq','camdt','camdq','ddt','ddq'},...
       {[1,-1; 1,1; -1,1; 1,1; 1,-1; 1,-1],...
        [-1,-1; 1,-1; 1,-1; -1,1; 1,1; -1,-1],...
        [1,-1; 1,1; -1,-1; -1,1; -1,-1; 1,-1],...
        [1,-1; -1,-1;-1,1; -1,-1; 1,1; -1,1],...
        [-1,-1; -1,-1; -1,-1; -1,1; -1,-1; 1,1],...
        [1,-1; 1,1; 1,-1; 1,-1; 1,1; 1,-1]});
  %
end
end
  sign_vec = sign_vec(vname);
  sign_vec = sign_vec(ii,:)
  marker = containers.Map(vnames,{'k-','k-.','k-','k-.'});
  if ~exist('ifhold','var')
    ifhold = 0;
  end
  if ifhold
    hold on
  else
    figure('units','normalized','outerposition',[0 0 1 1]); 
  end 
%  casei='F_2000_SPCAM_m2005_3hrly_f09f09_branch'
  casei='F_2000_SPCAM_m2005_3hrly2'
%  diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' season]
%  dlat = 'lat2525';
  diri = ['/Users/g/archive/matlab/' casei '/atm/hist/' dlat '/' season]
  load([diri '/mcs_clusters_1.mat'],'latC','nCl','nt4Cl'); % lat zone: trop(0),midlat(1),polar(2)
  load([diri '/mcs_cluster_parm.mat'])
%  ylimrange = [0 hp]
  ipm = parm.nz; %maximum hybrid level
  hp = flipud(parm.lev);
  hp = hp(1:ipm);
  ilndocn = ilndocn + 1; % +1 for actual index
  lndocnstr = {'ocn','lnd','lndocn','nocond'};
  latzone = {'trop','midlat','pole','global'};
  diro = ['/Users/g/archive/matlab/' casei '/figure/' dlat '/' season '/' lndocnstr{ilndocn} '/' latzone{ilatzone+1} '/']
  for iv=1:numel(vname1)
%%{
    if iv==2
      vtmp1 = vtmp;
    end
%%}
    eval(sprintf('load(''%s/mcs_cluster_stats/mcs_stats_%s'',''%s'',''%s'');',diri,lndocnstr{ilndocn},vname1{iv},'iclndocn2'));
    if strcmp(latzone{ilatzone+1},'global')
      ilatC = 1:numel(iclndocn2); % entire cluster over land without lat conditions
    else
      ilatC = find(latC(iclndocn2)==ilatzone); %clusters in the tropics(0) midlat(1)
      % latC: lat zone condition for all clusters, iclndocn2: cluster index
      % that are over lnd/ocn, ilatzone: lat zone condition
    end
    if dostatdata % use lifetime data, otherwise use original data without LS ave
      eval(sprintf('vtmp = %s.mcscomposite;',vname1{iv}));
  %%{
      if iv==2
        vtmp.data = vtmp1.data-vtmp.data;
      end
  %%}
      dat = squeeze(vtmp.data(:,ilife,ilatC));
    else
      eval(sprintf('load(''%s/mcs_cluster_var/mcs_cluster_%s'');',diri,vname1{iv}));
      eval(sprintf('data = %s;',vname1{iv}));
  %  ilatC = iclndocn2 & latC==ilatzone
      dat = [];
      for ic=1:numel(iclndocn2)
        % concat all the samples of ill into a matrix of (1:nz,1:allsamp)
        if ilife==1
          dat = cat(2,data{iclndocn2(ic)}{1},dat);
        elseif ilife==10 
          dat = cat(2,data{iclndocn2(ic)}{round(nt4Cl(iclndocn2(ic))/2)},dat);
        elseif ilife==20 
          dat = cat(2,data{iclndocn2(ic)}{end},dat);
        else
          error('ilife index out of range')
        end
      end
      dat = squeeze(dat(:,ilatC));
    end
  end
  if ismember(stats,{'eof'})
    [svec, variance, pc, pcmean, varmean] = svdVar(dat'); % start of lifetime
    if dosevsigntest
      plot(pcmean(1)*svec(:,1)*scale(vname1{1}),'r'); hold on
      plot(pcmean(2)*svec(:,2)*scale(vname1{1}),'b'); pause
      legend('EOF1','EOF2')
    end
%{
    if ilife==1
      signflip = [ 1 1 1];
    elseif ilife==numel(timerange)
      signflip = [ 1 -1 -1];
    end
%}
    varitmp = variance/sum(variance);
%   [svec, variance] = svdVar(squeeze(vtmp.data(:,end,:))'); % end of lifetime
%    for i = 1:numel(imode)
%      vtmp1 = vtmp1 + reshape(scale(ivv(vname))*vtmp.pc_mean(imode(i))*vtmp.svec(:,imode(i)),ipm,maxnt);
%    end
    marker_eof={'k-','k-.','k:','k.'};
    for i = 1:imode
      svec(:,i) = svec(:,i)*sign_vec(i); % this will go to function output
      pc(:,i) = pc(:,i)*sign_vec(i);
%      plot(scale(vname)*svec(:,i)*mean(pcmean(i)+pc(:,i)),hp,marker_eof{i},'linewidth',lw); hold on % the mean can be switched to median
      plot(svec(1:ipm,i),hp,marker_eof{i},'linewidth',lw); hold on
      lg{i} = (['EOF' num2str(i) ' (%Var=' num2str(round(varitmp(i)*100)) ')']);
    end
    plot(zeros(1,ipm),hp,'k','linewidth',0.5); hold on
    ll = legend(lg{:});
    ll.Box='off';
%{
    [~,loads] = sign_flip({pc(:,1:imode),svec(:,1:imode)},squeeze(vtmp.data(:,ilife,:))');
    [rvec varnew] = rotate_pca(loads{2}(:,1:imode),loads{1}(:,1:imode)); % (svec, pc)
    figure;
    for i = 1:imode
      plot(rvec(:,i),hp,marker{i},'linewidth',5); hold on
      lgr{i} = (['%Var = ' num2str(round(varnew(i)/sum(varnew)*100))]);
    end
    disp(['fractional variance for rotated pca modes = ' num2str(varnew/sum(varnew)*100)])
    legend(lgr{:})
%}
    xlim([-1 1])
%    ylim(ylimrange) 
    xlabel('EOF')
%    ylabel('height [km]')
    ylabel('hybrid level','FontSize',16)
    set(gca,'YTick',linspace(0,1000,6));
    set(gca,'YTickLabel',linspace(0,1000,6));
    ylim([0 1000])
    set(gca,'Ydir','reverse');
    set(gca,'Fontsize',44)
pause(3)
  elseif ismember(stats,{'fcbt'}) % fast chebyshev transform (1st kind) for non-periodic function
    zq = cgl(hp,ipm); % Chebyshev-Gauss-Lobatto points
    vq = interp1(hp,squeeze(scale(vname1{iv})*vtmp.data(:,ilife,:)),zq); % solutions on Chebyshev-Gauss-Lobatto point 
%{ 
    %consistency check   
    vq = interp1(hp,cos(hp*pi/hp),zq);
    plot(zq,vq,'b*'); hold on; pause; plot(zq,cos(zq*pi/hp),'r*'); pause
%}
    for isamp = 1:size(vtmp.data,3)
%{
      plot(vq(:,isamp)); hold on; 
      plot(squeeze(vtmp.data(:,ilife,isamp)),'r'); pause
%}
      vn(:,isamp) = fcht(vq(:,isamp)); % the coefficients of Chebyshev polynomials
%      plot(vn);hold on; pause
    end
    if iv==1
      plot([0:ipm-1],median(vn,2),marker(vname1{iv}),'linewidth',lw); hold on
    else
      plot([0:ipm-1],median(vn,2),'r:','linewidth',lw); hold on
    end
    xlim([0 ipm-1])
    xlabel('degree of polynomial')
    ylabel('amplitude')
    set(gca,'Fontsize',44)
    set(gca,'Xtick',[0:ipm-1])

  elseif ismember(stats,{'fft'})
    nz = size(vtmp.data,1);
    for isamp = 1:size(vtmp.data,3)
      V = fft(scale(vname1{iv})*vtmp.data(:,ilife,isamp)); % V = FT of variable v
      Vmag(:,isamp) = abs([V(1); 2*V(2:floor(nz/2)); V(floor(nz/2)+1)]); % only half + 1 needed for real number fft, V(1) = sum(vtmp.data(:,ilife,isamp))
%{
plot(vtmp.data(:,ilife,isamp),'r');hold on
plot(Vmag(:,isamp),'k'); pause; hold off
%}
      % Fourier Transform of v, plot a histogram of each freq later
    end
    wavelength = [inf hp./(1:floor(ipm/2))];
    plot(median(Vmag,2));
    xticklabels(wavelength);
%    ylim([0 5e-3])
%    ylim([0 5e-3])
  end
  if dosave
    if ilife == 1
      fname = ['EOF_' vname '_0'];
    elseif ilife==10
      fname = ['EOF_' vname '_10'];
    elseif ilife==20
      fname = ['EOF_' vname '_end'];
    end
    saveas(gcf,[diro '/' fname '.png']);
    crop([diro '/' fname '.png']);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function zq = cgl(zint,n)
  % Purpose: calc Chebyshev-Gauss-Lobatto points on an arbitrary interval other than [-1 1]
  zq = fliplr((zint(end)+zint(1))*0.5 + (zint(end)-zint(1))*0.5*cos((2*[1:n]-1)*pi*0.5./n));
 
 function B=fcgltran(A,direction)
% Fast Chebyshev Transform 
%
% Performs the fast transform of data sampled at the
% Chebyshev-Gauss-Lobatto Nodes x=cos(k*pi/N);
%
% A - original data in columns
% B - transformed data in columns
% direction - set equal to 1 for nodal to spectral
%             anything else for spectral to nodal
%
% Written by Greg von Winckel 03/08/04  
% Contact: gregvw@chtm.unm.edu

  [N,M]=size(A);
  
  if direction==1 % Nodal-to-spectral
      F=ifft([A(1:N,:);A(N-1:-1:2,:)]);
      B=([F(1,:); 2*F(2:(N-1),:); F(N,:)]);
  else            % Spectral-to-nodal
      F=fft([A(1,:); [A(2:N,:);A(N-1:-1:2,:)]/2]);
      B=(F(1:N,:));
  end 


function a = fcht(v)
% https://www.mathworks.com/matlabcentral/fileexchange/44030-fast-chebyshev-transform
%FCHT    Fast Chebyshev transform
% FCHT(V) computes the Chebyshev transform of a N+1 by 1 array V.  If V
% corresponds to a function evaluated at the Chebyshev points
% cos(pi*(0:N)/N), then V is interpolated by a linear combinations of the
% Chebyshev polynomials with weights given by FCHT(V).
% 
% 
% Example:
% Approximate f(x) = exp(x) over [-1,1] as a linear combination of the
% first three Chebyshev polynomials.
% 
% x = cos(pi*(0:2)/2); % establish 3 Chebyshev grid points
% 
% V = exp(x); % evaluate f(x) at Chebyshev grid points
% 
% a = fcht(V);
% xx = linspace(-1,1); % create dense grid over domain
% g = a(1)*1 + a(2)*xx + a(3)*(2*xx.^2 - 1); % sum the first three Chebyshev
%               % polynomials with respect to their corresponding weights
% plot(xx,exp(xx),xx,g); % visualize the approximation

v = v(:);
N = length(v) - 1;
v = [v; flipud(v(2:N))];

a = real(fft(v))/N;
a = [a(1)/2; a(2:N); a(N+1)/2];








%{


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%
casei = 'F_2000_SPCAM_m2005_3hrly2';
%season = 'DJF';
season = 'JJA';
lndocn = {'ocn','lnd'};
%vname = 'spdt';
%vname = 'spdq';
%vname = 'camdt';
vname = 'camdq';

eval(sprintf(['svec_' vname '= svec';])) 
eval(sprintf(['pc_' vname  '= pc';])) % do this for camdt after rerun
% {'ocn 0%','ocn 100%','lnd 0%','lnd 100%'}
% svec_spdt = svec; pc_spdt = pc; % do this for camdt after rerun
% svec_camdt = svec; pc_camdt = pc; % do this for camdt after rerun
ii=1;
nmode=2; % EOF3 variance is approx 10% or less for both camdt/camdq or spdt/spdq
for ilnd = [0 1]
  for ilife = [1 10 20]
%    [R{ii} P{ii}] = corr(pc_spdt{ii}(:,1:nmode),pc_camdt{ii}(:,1:nmode));
    [R{ii} P{ii}] = corr(pc_spdq{ii}(:,1:nmode),pc_camdq{ii}(:,1:nmode));
    ii=ii+1;
end; end
diro = ['/Users/g/archive/matlab/F_2000_SPCAM_3hrly2/figure/' dlat '/' season '/' latzone{ilatzone+1} '/'];
%save([diro 'corr_spdt_camdt_pc.mat'],'R','P','pc_spdt','pc_camdt');
save([diro '/corr_spdq_camdq_pc.mat'],'R','P','pc_spdq','pc_camdq');
% R =
%
%    0.1242    0.1436    0.0539
%    0.4438   -0.4749    0.1007
%   -0.3468   -0.2416    0.1416
%
%
% P =
%
%    0.0000    0.0000    0.0016
%    0.0000    0.0000    0.0000
%    0.0000    0.0000    0.0000
%
% [R P] = corr(pc_spdt{2}(:,1:3),pc_camdt{2}(:,1:3))
% R =
%
%    0.5026    0.1573    0.1607
%    0.4484   -0.5318    0.1844
%    0.2418    0.2475   -0.2369
%
%
% P =
%
%   1.0e-19 *
%
%    0.0000    0.2257    0.0350
%    0.0000    0.0000    0.0000
%    0.0000    0.0000    0.0000
%
[R P] = corr(pc_spdt{3}(:,1:3),pc_camdt{3}(:,1:3))

R =

    0.5026    0.1573    0.1607
    0.4484   -0.5318    0.1844
   -0.2418   -0.2475    0.2369


P =

   1.0e-19 *

    0.0000    0.2257    0.0350
    0.0000    0.0000    0.0000
    0.0000    0.0000    0.0000

[R P] = corr(pc_spdt{4}(:,1:3),pc_camdt{4}(:,1:3))

R =

    0.1589    0.1375    0.0571
    0.4728   -0.3537    0.0154
    0.0192   -0.4696    0.1468


P =

    0.0000    0.0000    0.0327
    0.0000    0.0000    0.5655
    0.4724    0.0000    0.0000
% [R P] = corr(pc_spdt{5}(:,1:3),pc_camdt{3}(:,1:3))
% R =
%
%    0.1589    0.1375    0.0571
%    0.4728   -0.3537    0.0154
%   -0.0192    0.4696   -0.1468
%
% P =
%
%    0.0000    0.0000    0.0327
%    0.0000    0.0000    0.5655
%    0.4724    0.0000    0.0000
% [R P] = corr(pc_spdt{6}(:,1:3),pc_camdt{4}(:,1:3))
% R =
%
%    0.3275    0.2374    0.4340
%    0.6010   -0.3692    0.1385
%    0.2036    0.4680   -0.2432
%
%
% P =
%
%   1.0e-06 *
%
%    0.0000    0.0000    0.0000
%    0.0000    0.0000    0.2035
%    0.0000    0.0000    0.0000
%
% [R P] = corr(pc_spdq{1}(:,1:3),pc_camdq{1}(:,1:3))
% R =
%
%   -0.0891   -0.0019   -0.1490
%   -0.3344    0.0338   -0.0738
%    0.2377   -0.0993    0.0698
%
% P =
%
%    0.0000    0.9101    0.0000
%    0.0000    0.0480    0.0000
%    0.0000    0.0000    0.0000
% [R P] = corr(pc_spdq{2}(:,1:3),pc_camdq{2}(:,1:3))
% R =
%
%    0.0668   -0.1489   -0.2729
%   -0.3956    0.0161   -0.0722
%    0.1631   -0.1530   -0.0573
%
% P =
%
%    0.0001    0.0000    0.0000
%    0.0000    0.3481    0.0000
%    0.0000    0.0000    0.0008
%
% [R P] = corr(pc_spdq{5}(:,1:3),pc_camdq{3}(:,1:3))
% R =
%
%    0.0454   -0.0540   -0.1624
%   -0.4947   -0.1268   -0.1082
%    0.0558    0.1482    0.0726
%
% P =
%
%    0.0897    0.0437    0.0000
%    0.0000    0.0000    0.0001
%    0.0372    0.0000    0.0066
%
% [R P] = corr(pc_spdq{6}(:,1:3),pc_camdq{4}(:,1:3))
% R =
%
%    0.1025   -0.2304   -0.1517
%   -0.5753   -0.0416   -0.0567
%    0.2028    0.0796    0.1710
%
% P =
%
%    0.0001    0.0000    0.0000
%    0.0000    0.1199    0.0341
%    0.0000    0.0029    0.0000
%%%%%%%%%%%%% Discrete Chebyshev Transform %%%%%%%%%%%%%%%%%%%%%%
tq='dt'; ipm=20; season = 'JJA'; 
ii=1;
for ilnd= [0 1]
  for ilife = [1 20]
    coef_sp{ii} = plot_mcs_profile({['sp' tq]},'fcbt',ilnd,season,0,ilife); hold on;
    coef_cam{ii} = plot_mcs_profile({['cam' tq]},'fcbt',ilnd,season,1,ilife); 
    coef_sp_cam{ii} = plot_mcs_profile({['sp' tq],['cam' tq]},'fcbt',ilnd,season,1,ilife);
    ii=ii+1;
    switch tq
    case 'dq' 
      ylim([-15 15]);
    case 'dt' 
      ylim([-15 15]);
    end
    if ilife == 1 & ilnd == 1 
      hl = legend('SPCAM','CAM','SPCAM - CAM');
      set(hl,'box','off')
    end
    plot([0:ipm-1],zeros(1,ipm),'k');
    lndocn = {'ocn','lnd'};
    if ilife == 1
      fname = ['FCT_' tq '_0'];
    elseif ilife==20
      fname = ['FCT_' tq '_end'];
    end
   diro=['/Users/g/archive/matlab/F_2000_SPCAM_m2005_3hrly_f09f09_branch/figure/'];
    diro=[diro season '/' lndocn{ilnd+1} '/'];
    saveas(gcf,[diro '/' fname '.png']);
    crop([diro '/' fname '.png']);
  end
end



%}
