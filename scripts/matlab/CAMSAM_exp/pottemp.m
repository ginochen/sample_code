function th = pottemp(T,p)
   P0   = 100000; % ref pressure for potential temp
   Rd   = 287;
   Cp   = 1004; %. or 1005.7 % specific heat dry air [J/kg/K]
   RCP  = Rd/Cp; 
   th     = T.*(P0./p).^RCP; % potential temp of the parcel 
