function res = dvdt_narrow(u,t,b,c11,c12,c22,F,d1,B,p)
% Approximate v_t = b*((c11*u_x+c12*u_y)_x+(c12*u_x+c22*u_y)_y)+F using
% 2p-th order narrow stencil d1 & B

% (c11*u_x)_x
Dxx_u = zeros(size(u));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        Dxx_u = Dxx_u + B(m+p+1,n+p+1)*circshift(c11,-m,2).*circshift(u,-n,2);
    end
end

% u_y
Dy_u = zeros(size(u));
for n = -p:p
    Dy_u = Dy_u + d1(n+p+1)*circshift(u,-n,1);
end
% c12*u_y
Dy_u = c12.*Dy_u;
% (c12*u_y)_x
DxDy_u = zeros(size(u));
for n = -p:p
    DxDy_u = DxDy_u + d1(n+p+1)*circshift(Dy_u,-n,2);
end

% u_x
Dx_u = zeros(size(u));
for n = -p:p
    Dx_u = Dx_u + d1(n+p+1)*circshift(u,-n,2);
end
% c12*u_x
Dx_u = c12.*Dx_u;
% (c12*u_x)_y
DyDx_u = zeros(size(u));
for n = -p:p
    DyDx_u = DyDx_u + d1(n+p+1)*circshift(Dx_u,-n,1);
end

% (c22*u_y)_y
Dyy_u = zeros(size(u));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        Dyy_u = Dyy_u + B(m+p+1,n+p+1)*circshift(c22,-m,1).*circshift(u,-n,1);
    end
end

% b*((c11*u_x+c12*u_y)_x+(c12*u_x+c22*u_y)_y)+F
res = b.*(Dxx_u+DxDy_u+DyDx_u+Dyy_u)+F(t);

end
