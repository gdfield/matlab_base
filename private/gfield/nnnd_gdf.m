function d = nnnd_gdf(com1, com2, angle1, angle2, sigmas1, sigmas2)
%
% usage:  d = nnnd_gdf(com1, com_2, angle1, angle2, sigmas1, sigmas2)
%
% arguments:    com1 - x-y vector, center of mass of cell 1
%               angle1 - angle from fit to the sta            
%               sigmas1 = 2d vector, sigmax and sigmay for cell 1
%               com2 - x-y vector, center of mass of cell 2
%               angle2 - angle from fit to the sta            
%               sigmas2 = 2d vector, sigmax and sigmay for cell 2
%
% outputs:     d - nnnd between two cells
%
%
% 2018-04 GDF



% assumes format = [xctr yctr theta sd_max sd_min] 
v = com1-com2;

t = atan2(v(2),v(1));

t0 = angle1;
%sd = a(4:5);
sa = norm([sigmas1(1)*cos(t-t0) sigmas1(2)*sin(t-t0)]);
        
t0 = angle2;
%sd = b(4:5);        
sb = norm([sigmas2(1)*cos(t-t0) sigmas2(2)*sin(t-t0)]);

sm = (sa+sb);
d = norm(v)/sm;