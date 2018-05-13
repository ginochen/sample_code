% Xk_orig = Xk_anom + Xk_mean 
load ~/archive/matlab/var/var/Q_ZMincUW_UWexcZM_basinwise-svd.mat
docam=1;
if docam
  pc = pc_z; 
  pcmean = pcmean_z; 
  svec = svec_z;
  var = var_z;
else % do spcam
  pc = pc_spz; 
  pcmean = pcmean_spz; 
  svec = svec_spz;
  var = var_spz;
end
%nmode = size(svec{ib},2); % all modes
nmode = 3;
nz = size(svec{1},1);
dz = 17;
zvec = [nz-dz:nz]; % use the surface (nz) to dz+1-th (nz-dz) pressure level, notice pressure increases with index
xmin = -10;
xmax = 10;
nbin = 200;
ii=1;
for ib=1:4;
  ntsamp = size(pc{ib},1);
  fracvar = var{ib}/sum(var{ib});
  for im=1:nmode
    xkb = pcmean{ib}(im)*svec{ib}(zvec,im);
    for iz = 1:dz+1
      xk_orig{ib}(iz,:) = xkb(iz) + pc{ib}(:,im)*svec{ib}(zvec(iz),im);  
    end
    subplot(4,nmode,ii)
    [vprob, vbincenters] = binVarWithHeight(xk_orig{ib},xmax,xmin,nbin,dz+1);
    contourf(vbincenters,1:dz+1,flipud(vprob),1000,'linestyle','none')
    hold on; plot(flipud(xkb),1:dz+1,'r-')
    text(vbincenters(end-2),2,num2str(fracvar(im)))
    colormap(cmap(5))
    display(im)
    %display(['fractional variance = ' num2str(fracvar(im))])
    ii=ii+1;
  end
end

