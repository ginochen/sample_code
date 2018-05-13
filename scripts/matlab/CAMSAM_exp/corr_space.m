
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% space-correlated, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
iilev=[1:30];
rho.iilev = iilev;
[rho.e_precc     p.e_precc]         = corrcoef_2d_vert_nzero(var{31},var{4},ilat,ilon,itimes,iilev);
[rho.eabs_precc  p.ebs_precc]       = corrcoef_2d_vert_nzero(var{32},var{4},ilat,ilon,itimes,iilev);
[rho.e_preccdzm  p.e_preccdzm]      = corrcoef_2d_vert_nzero(var{31},var{10},ilat,ilon,itimes,iilev);
[rho.eabs_preccdzm  p.ebs_preccdzm] = corrcoef_2d_vert_nzero(var{32},var{10},ilat,ilon,itimes,iilev);
[rho.e_precsh  p.e_precsh]          = corrcoef_2d_vert_nzero(var{31},var{11},ilat,ilon,itimes,iilev);
[rho.eabs_precsh  p.eabs_precsh]     = corrcoef_2d_vert_nzero(var{32},var{11},ilat,ilon,itimes,iilev);
[rho.e_relhum    p.e_relhum]        = corrcoef_2d_vert(var{31},var{18},ilat,ilon,itimes,iilev);
[rho.eabs_relhum p.eabs_relhum]     = corrcoef_2d_vert(var{32},var{18},ilat,ilon,itimes,iilev);
[rho.e_omegat     p.e_omegat]       = corrcoef_2d_vert(var{31},var{15},ilat,ilon,itimes,iilev);
[rho.eabs_omegat  p.eabs_omegat]     = corrcoef_2d_vert(var{32},var{15},ilat,ilon,itimes,iilev);
[rho.e_cape     p.e_cape]           = corrcoef_2d_vert_nnan(var{31},var{8},ilat,ilon,itimes,iilev);
[rho.eabs_cape  p.eabs_cape]         = corrcoef_2d_vert_nnan(var{32},var{8},ilat,ilon,itimes,iilev);
[rho.e_cin     p.e_cin]             = corrcoef_2d_vert_nnan(var{31},var{9},ilat,ilon,itimes,iilev);
[rho.eabs_cin  p.eabs_cin]           = corrcoef_2d_vert_nnan(var{32},var{9},ilat,ilon,itimes,iilev);
savefig(gcf,'corr_e_cape')
