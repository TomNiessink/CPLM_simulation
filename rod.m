function [H,angH] = rod(r,l,theta,Gridsize_m,Gridsize_px) 
% Function which generates a rod-shaped crystal with
% Length l (m), radius r (m), angle theta to the optical axis (rad),
% and a defined grid (Gridsize_m, m and Gridsize_px, px)

H = zeros (Gridsize_px); 
res = Gridsize_m/Gridsize_px;

crys = zeros(round(2*r/res)+1,round(l/res));
angles = zeros(round(2*r/res)+1,round(l/res));
angles(:,:) = theta;

x=-r:res:r;
h_x=2*r*sin(acos(x/r));
for i=1:round(l/res)
crys(:,i)=h_x;
end
crys = imrotate(crys,rad2deg(theta),'bilinear');
angles = imrotate(angles,rad2deg(theta),'bilinear');

szcrys = size(crys);
szH = size(H);
angH=H;

H((szH(1)/2-szcrys(1)/2):(szH(1)/2-szcrys(1)/2+szcrys(1)-1),...
    (szH(2)/2-szcrys(2)/2):(szH(2)/2-szcrys(2)/2+szcrys(2))-1)=crys;

angH((szH(1)/2-szcrys(1)/2):(szH(1)/2-szcrys(1)/2+szcrys(1)-1),...
    (szH(2)/2-szcrys(2)/2):(szH(2)/2-szcrys(2)/2+szcrys(2))-1)=angles;



