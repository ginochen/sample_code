function [ WRF_DBZ, WRF_DBZ_MAX, Z_E ] = wrf_dbz(PRES,TK,QRAIN,QGRAUP,QSNOW,QVAPOR)
% https://www.ncl.ucar.edu/Support/talk_archives/2012/att-0941/wrf_user_dbz.f
%     This routine computes equivalent reflectivity factor (in dBZ) at
%     each model grid point.  In calculating Ze, the RIP algorithm makes
%     assumptions consistent with those made in an early version
%     (ca. 1996) of the bulk mixed-phase microphysical scheme in the MM5
%     model (i.e., the scheme known as "Resiner-2").  For each species:
%
%     1. Particles are assumed to be spheres of constant density.  The
%     densities of rain drops, snow particles, and graupel particles are
%     taken to be rho_r = rho_l = 1000 kg m^-3, rho_s = 100 kg m^-3, and
%     rho_g = 400 kg m^-3, respectively. (l refers to the density of
%     liquid water.)
%
%     2. The size distribution (in terms of the actual diameter of the
%     particles, rather than the melted diameter or the equivalent solid
%     ice sphere diameter) is assumed to follow an exponential
%     distribution of the form N(D) = N_0 * exp( lambda*D ).
%
%     3. If ivarint=0, the intercept parameters are assumed constant
%     (as in early Reisner-2), with values of 8x10^6, 2x10^7, 
%    and 4x10^6 m^-4, for rain, snow, and graupel, respectively.
%    If ivarint=1, variable intercept parameters are used, as 
%    calculated in Thompson, Rasmussen, and Manning (2004, Monthly
%    Weather Review, Vol. 132, No. 2, pp. 519-542.)
%
%     4. If iliqskin=1, frozen particles that are at a temperature above
%     freezing are assumed to scatter as a liquid particle.
%
%     More information on the derivation of simulated reflectivity in
%     RIP can be found in Stoelinga (2005, unpublished write-up).
%     Contact Mark Stoelinga (stoeling@atmos.washington.edu) for a copy.
%
% https://sourceforge.net/p/vapor/mailman/message/27509953/
% Python program to calculate radar reflectivity using WRF variables
% Copied from NCL/Fortran source code wrf_user_dbz.f
% Based on work by Mark Stoellinga, U. of Washington
% Inputs:  P, PB, QRAIN, QGRAUP, QSNOW, T, QVAPOR
% Outputs: WRF_DBZ (3d) and WRF_DBZ_MAX (2d)
% WRF_DBZ_MAX is a 2D variable, the maximum over vertical columns of WRF_DBZ
% parameters ivarint and iliqskin defined as in original code, either 0 or 1.
% If iliqskin=1, frozen particles above freezing are assumed to scatter as a liquid particle
% If ivarint=0, the intercept parameters are assumed constant with values of 
%  8*10^6, 2*10^7, and 4*10^6, for rain, snow, and graupel respectively.
% If ivarint=1, intercept parameters are used as calculated in Thompson, Rasmussen, and Manning,
% 2004 Monthly Weather Review, Vol 132, No. 2, pp.519-542
%  By default they are given value 0, but they can be changed by editing the following 2 lines:
iliqskin    = 1; % default is zero
ivarint     = 0;
%c           = 2.0/7.0;

% Constants used to calculate variable intercepts
R1          = 1.e-15;
RON         = 8.e6;
RON2        = 1.e10;
SON         = 2.e7;
GON         = 5.e7; % there's a typo in python code as 5.37
RON_MIN     = 8.e6;
RON_QR0     = 0.0001;
RON_DELQR0  = 0.25*RON_QR0;
RON_CONST1R = (RON2-RON_MIN)*0.5;
RON_CONST2R = (RON2+RON_MIN)*0.5;

% Constant intercepts
RN0_R       = 8.e6;
RN0_S       = 2.e7;
RN0_G       = 4.e6;

% Other constants
GAMMA_SEVEN = 720.;
RHOWAT      = 1000.;
RHO_R       = RHOWAT;
RHO_S       = 100.;
RHO_G       = 400.;
ALPHA       = 0.224;
CELKEL      = 273.15;
PI          = 3.141592653589793;
RD          = 287.04;

 
% #calculate Temp. in Kelvin
% PRES = P+PB
% TK = (T+300.)*numpy.power(PRES*.00001,c)
 
% Force Q arrays to be nonnegative:
QVAPOR = max(QVAPOR,0.0);
QSNOW  = max(QSNOW,0.0);
QRAIN  = max(QRAIN,0.0);
if (isempty(QGRAUP)~=0)
   QGRAUP = max(QGRAUP,0.0); % don't know why they use QGRAUP = 0.0 here
else
   QGRAUP = 0.0;
end 
%TestT = less(TK,CELKEL); % find the index that are below freezing
%QSNOW = where(TestT,QRAIN, QSNOW); % if below freezing use QRAIN, else use QSNOW
%QRAIN = where(TestT,0.0, QRAIN);
iless = find(TK<CELKEL); % below freezing
imore = find(TK>=CELKEL); % above freezing
QSNOW(iless) = QRAIN(iless);
QRAIN(iless) = 0;
 

VIRTUAL_T = TK.*(0.622 + QVAPOR)./(0.622*(1.+QVAPOR));
RHOAIR = PRES./(RD*VIRTUAL_T);
 
FACTOR_R    = GAMMA_SEVEN*1.e18*power((1./(PI*RHO_R)),1.75); % const
FACTOR_S    = GAMMA_SEVEN*1.e18*power((1./(PI*RHO_S)),1.75)*power((RHO_S/RHOWAT),2.*ALPHA); % const
FACTOR_G    = GAMMA_SEVEN*1.e18*power((1./(PI*RHO_G)),1.75)*power((RHO_G/RHOWAT),2.*ALPHA);

% Adjust factor for brightband, where snow or graupel particle
%      scatters like liquid water (alpha=1.0) because it is assumed to
%      have a liquid skin. 
if (iliqskin == 1);
   %FACTORB_S = float32(where(TestT,FACTOR_S,FACTOR_S/ALPHA)); % substitute TRUE with FACTOR_S, FALSE with FACTOR_S/ALPHA
   %FACTORB_G = float32(where(TestT,FACTOR_G, FACTOR_G/ALPHA));
   FACTORB_S(iless) = FACTOR_S;
   FACTORB_G(iless) = FACTOR_G;
   FACTORB_S(imore) = FACTOR_S/ALPHA;
   FACTORB_G(imore) = FACTOR_G/ALPHA;
else
   FACTORB_S = FACTOR_S;
   FACTORB_G = FACTOR_G;
end

% Calculate variable intercept parameters
if (ivarint == 1);
   TEMP_C = min(-.001, TK-CELKEL); 
   SONV   = min(2.e8, 2.e6*exp(-.12*TEMP_C));
   imore  = find(QGRAUP>R1);
   iless  = find(QGRAUP<=R1);
   GONV(imore) = 2.38*power(PI*RHO_G/(RHOAIR(imore).*QGRAUP(imore)),0.92);
   GONV(iless) = GON(iless);
   GONV(imore) = max(1.e4,min(GONV(imore),GON(imore)));
   imore = find(QRAIN>R1);
   iless = find(QRAIN<=R1);
   RONV(imore) = RON_CONST1R*tanh((RON_QR0-QRAIN(imore))/RON_DELQR0) + RON_CONST2R;
   RONV(iless) = RON2; 
   %TestG = less(R1, QGRAUP);
   %GONV = where(TestG,2.38*power(PI*RHO_G/(RHOAIR*QGRAUP),0.92),GON);
   %GONV = where(TestG,maximum(1.e4,minimum(GONV,GON)),GONV);
   %RONV = where(greater(QRAIN,R1),RON_CONST1R*tanh((RON_QR0-QRAIN)/RON_DELQR0) + RON_CONST2R,RON2);
else;
   GONV = RN0_G;
   RONV = RN0_R;
   SONV = RN0_S;
end

% Total equivalent reflectivity factor (Z_E, in mm^6 m^-3) is
%      the sum of z_e for each hydrometeor species: 
Z_E = FACTOR_R*power(RHOAIR.*QRAIN,1.75)./power(RONV,0.75) + ...
      FACTORB_S'.*power(RHOAIR.*QSNOW,1.75)./power(SONV,0.75) + ...
      FACTORB_G'.*power(RHOAIR.*QGRAUP,1.75)./power(GONV,0.75);

% Adjust small values of Z_E so that dBZ is no lower than -30
Z_E = max(Z_E, 0.001);

WRF_DBZ = 10.0*log10(Z_E); % convert Z_E to DBZ
WRF_DBZ_MAX = max(WRF_DBZ); % check if the max dimension is 1
%WRF_DBZ_MAX = amax(WRF_DBZ,axis=0); % amax := array max 
