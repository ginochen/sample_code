function pr = rotate_storm(uvstorm)
% Rotate the initial point (pi) 
%pi=[-3,-5]; 
if uvstorm(1)<0; signs=-1; else; signs=1; end; 
[th] = degrees(uvstorm',[0 1]'); % rotate degrees wrt y-coord [0 1]
for ix = 
  for iy = 
pr = rotate_coord(signs*th,[3,5])
