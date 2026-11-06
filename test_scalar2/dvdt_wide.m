function res = dvdt_wide(u,t,b,c11,c12,c22,F,d1,p)
% Approximate v_t = b*((c11*u_x+c12*u_y)_x+(c12*u_x+c22*u_y)_y)+F using
% 2p-th order narrow stencil d1 & wide stencil d1^2

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
% c11*u_x+c12*u_y
sum1 = c11.*Dx_u+c12.*Dy_u;
% (c11*u_x+c12*u_y)_x
for n = -p:p
    res = res + d1(n+p+1)*circshift(sum1,-n,2);
end
% c12*u_x+c22*u_y
sum2 = c12.*Dx_u+c22.*Dy_u;
% (c21*u_x+c22*u_y)_y
for n = -p:p
    res = res + d1(n+p+1)*circshift(sum2,-n,1);
end

% b*((c11*u_x+c12*u_y)_x+(c12*u_x+c22*u_y)_y)+F
res = b.*res+F(t);

end
