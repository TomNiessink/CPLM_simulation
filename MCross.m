function [H,angH] = MCross(r,ne,no,Gridsize_m,Gridsize_px) 




res = Gridsize_m/Gridsize_px;
lb = -1*Gridsize_m/2;
rb = Gridsize_m/2;

x = lb:res:rb; %m
y = lb:res:rb; %m
z = lb:res:rb; %m

r_ij = zeros(size(x,2),size(y,2));
d_ij = zeros(size(x,2),size(y,2));
ang_ij = zeros(size(x,2),size(y,2));
ne_p_ijk = zeros(size(x,2),size(y,2),size(z,2));
dne_p_ijk = zeros(size(x,2),size(y,2),size(z,2));
r_ijk = zeros(size(x,2),size(y,2),size(z,2));
theta_ijk = zeros(size(x,2),size(y,2),size(z,2));

for i = 1:size(x,2)
    for j = 1:size(x,2)
        ang_ij(i,j) = atan(y(j)/x(i));
        for k = 1:size(z,2)
        r_ij(i,j) = sqrt(x(i)^2+y(j)^2);
        
        d_ij(i,j) = 2*sqrt(abs(r^2-r_ij(i,j)^2));
        r_ijk(i,j,k) = sqrt(z(k)^2+r_ij(i,j)^2);
    
       
        if r_ijk(i,j,k) > r
            ne_p_ijk(i,j,k)=0;
        else
            theta_ijk(i,j,k) = pi/2-atan((d_ij(i,j)/2-z(k))/r_ij(i,j));
            ne_p_ijk(i,j,k) = 1/sqrt((cos(theta_ijk(i,j,k))^2/no^2)+(sin(theta_ijk(i,j,k))^2/ne^2));
            dne_p_ijk(i,j,k) = no-ne_p_ijk(i,j,k);
        end
        end
    end
end




H = sum(dne_p_ijk,3)*res;
H(isnan(H))=0;
angH = ang_ij;
angH(isnan(angH))=0;

surf(angH)
