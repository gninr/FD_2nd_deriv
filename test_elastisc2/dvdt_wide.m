function res = dvdt_wide(u,t,b,c,F,d1,p)
% Approximate (v1,v2)_t = b*((c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y)_x +
% (c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y)_y) + F using 2p-th order
% narrow stencil d1 & wide stencil d1^2

% u_x
Dx_u = zeros(size(u));
for n = -p:p
    Dx_u = Dx_u + d1(n+p+1)*circshift(u,-n,2);
end
% u_y
Dy_u = zeros(size(u));
for n = -p:p
    Dy_u = Dy_u + d1(n+p+1)*circshift(u,-n,1);
end

res = zeros(size(u));
% c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y
sum1 = c{1,1,1}.*Dx_u(:,:,1)+c{1,2,1}.*Dy_u(:,:,1)+c{1,1,2}.*Dx_u(:,:,2)+c{1,2,2}.*Dy_u(:,:,2);
% (c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y)_x
for n = -p:p
    res = res + d1(n+p+1)*circshift(sum1,-n,2);
end
% c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y
sum2 = c{2,1,1}.*Dx_u(:,:,1)+c{2,2,1}.*Dy_u(:,:,1)+c{2,1,2}.*Dx_u(:,:,2)+c{2,2,2}.*Dy_u(:,:,2);
% (c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y)_y
for n = -p:p
    res = res + d1(n+p+1)*circshift(sum2,-n,1);
end

% b.*((c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y)_x+(c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y)_y)+F
res = b.*res+F(t);

end
