load('var_PC1_11.mat','parm')
load var_PC1_MCSt.mat
idv = getVarIndex({'varCRM2D','varCAM1D'},parm);
docam=1;
if docam
  [svec var pc pcmean] = svdVar(squeeze((varCAM1D_mcsall(:,idv.varCAM1D.u,:))),1,1,0,1);
xmin = [-20 -20 -20];
xmax = [20 20 20];
ii=1;
else % do spcam
  [svec var pc pcmean] = svdVar(squeeze(mean(varCRM2D_mcsall(:,:,idv.varCRM2D.u,:),1)),1,1,0,1);
xmin = [-20 -20 -20];
xmax = [20 20 20];
ii=4;
end
%nmode = size(svec{ib},2); % all modes
nmode = 3;
nz = size(svec,1);
zvec = [1:nz]; % use the surface (nz) to dz+1-th (nz-dz) pressure level, notice pressure increases with index
zint = parm.zint/1000;
nbin = 200;
  ntsamp = size(pc,1);
  fracvar = var/sum(var);
  for im=1:nmode
    xkb = pcmean(im)*svec(zvec,im);
    for iz = 1:nz
      xk_orig(iz,:) = xkb(iz) + pc(:,im)*svec(zvec(iz),im);  
    end
    subplot(2,nmode,ii)
    [vprob, vbincenters] = binVarWithHeight(xk_orig,xmax(im),xmin(im),nbin,nz);
    contourf(vbincenters,zint,vprob,1000,'linestyle','none')
    hold on; plot(xkb,zint,'r-')
    text(vbincenters(end-2),2,num2str(fracvar(im)))
    colormap(cmap(5))
    display(im)
    %display(['fractional variance = ' num2str(fracvar(im))])
    ii=ii+1;
  end

