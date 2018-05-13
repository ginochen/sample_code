% Purpose: Assemble the ke spectrum for all selected timesteps
for it=11:10:500
   load('/projects/rsmas/kirtman/gchen/archive/matlab/figure/var/ke_spectrum/PC1_fft-basinwise' num2str(it) '.mat');
   for ib=1:4
   FFT{ib}(:,:,:,1:size(P1{ib},4)) = P1{ib} % P1(nx/2+1 (fft modes), nz (17 levels surface to 200mb), nlag+1 (7 lags), nlatlon (at the selected timestep)) 
             % iillt_zm_PC is the same for all timestep files
   
end
