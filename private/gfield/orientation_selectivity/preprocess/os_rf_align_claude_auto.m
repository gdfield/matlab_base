function A = os_rf_align_claude_auto(z, pol, row_c, col_c, theta, hw, scale)
% os_rf_align_claude_auto  Signed RF rotated so the preferred drift axis lies at 0 deg.
% z: robust-z spatial RF (rows x cols, row 1 = top of display). pol: +1 ON / -1 OFF
% (0/NaN treated as +1). (row_c, col_c): RF center in array coordinates. theta:
% preferred drift axis (deg, counterclockwise from rightward, y up = screen up;
% grating and STA coordinates assumed to match, GDF 2026-09-27). hw: half width
% (default 12 -> 25x25). Output A(v,u): row 1 = +hw (up), column 1 = -hw; bilinear
% interpolation, 0 outside the image; unit Frobenius norm. scale (default 1): stixels
% of the source image per output pixel (used to match RF size to the template).
% 2026-09-28 GDF + Claude
if nargin < 6 || isempty(hw), hw = 12; end
if nargin < 7 || isempty(scale), scale = 1; end
if isnan(pol) || pol == 0, pol = 1; end
[U, V] = meshgrid(-hw:hw, hw:-1:-hw);
c = cosd(theta); s = sind(theta);
col = col_c + scale * (U * c - V * s);
row = row_c - scale * (U * s + V * c);
A = interp2(pol * z, col, row, 'linear', 0);
nA = norm(A(:)); if nA > 0, A = A / nA; end
end
