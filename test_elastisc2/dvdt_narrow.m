function res = dvdt_narrow(u,t,b,c,F,d1,B,p)
% Approximate (v1,v2)_t = b*((c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y)_x +
% (c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y)_y) + F using 2p-th order
% narrow stencil d1 & B

% (c111*u1_x)_x & (c112*u2_x)_x
Dxx_u1 = zeros(size(u));
Dxx_u2 = zeros(size(u));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        Dxx_u1 = Dxx_u1 + B(m+p+1,n+p+1)*circshift(c{1,1,1},-m,2).*circshift(u(:,:,1),-n,2);
        Dxx_u2 = Dxx_u2 + B(m+p+1,n+p+1)*circshift(c{1,1,2},-m,2).*circshift(u(:,:,2),-n,2);
    end
end

% u_y
Dy_u = zeros(size(u));
for n = -p:p
    Dy_u = Dy_u + d1(n+p+1)*circshift(u,-n,1);
end
% c121*u1_y
Dy_u1 = c{1,2,1}.*Dy_u(:,:,1);
% c122*u2_y
Dy_u2 = c{1,2,2}.*Dy_u(:,:,2);
% (c121*u1_y)_x & (c122*u2_y)_x
DxDy_u1 = zeros(size(u));
DxDy_u2 = zeros(size(u));
for n = -p:p
    DxDy_u1 = DxDy_u1 + d1(n+p+1)*circshift(Dy_u1,-n,2);
    DxDy_u2 = DxDy_u2 + d1(n+p+1)*circshift(Dy_u2,-n,2);
end

% u_x
Dx_u = zeros(size(u));
for n = -p:p
    Dx_u = Dx_u + d1(n+p+1)*circshift(u,-n,2);
end
% c211*u1_x
Dx_u1 = c{2,1,1}.*Dx_u(:,:,1);
% c212*u2_x
Dx_u2 = c{2,1,2}.*Dx_u(:,:,2);
% (c211*u1_x)_y & (c212*u2_x)_y
DyDx_u1 = zeros(size(u));
DyDx_u2 = zeros(size(u));
for n = -p:p
    DyDx_u1 = DyDx_u1 + d1(n+p+1)*circshift(Dx_u1,-n,1);
    DyDx_u2 = DyDx_u2 + d1(n+p+1)*circshift(Dx_u2,-n,1);
end

% (c221*u1_y)_y & (c222*u2_y)_y
Dyy_u1 = zeros(size(u));
Dyy_u2 = zeros(size(u));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        Dyy_u1 = Dyy_u1 + B(m+p+1,n+p+1)*circshift(c{2,2,1},-m,1).*circshift(u(:,:,1),-n,1);
        Dyy_u2 = Dyy_u2 + B(m+p+1,n+p+1)*circshift(c{2,2,2},-m,1).*circshift(u(:,:,2),-n,1);
    end
end

% b.*((c111*u1_x+c121*u1_y+c112*u2_x+c122*u2_y)_x+(c211*u1_x+c221*u1_y+c212*u2_x+c222*u2_y)_y)+F
res = b.*(Dxx_u1+DxDy_u1+DyDx_u1+Dyy_u1+Dxx_u2+DxDy_u2+DyDx_u2+Dyy_u2)+F(t);

end
