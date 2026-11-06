function res = dvdt_wide(u,t,b,c,F,d1,p)
% Approximate v_t = b*(c*u_x)_x+F using 2p-th order wide stencil d1^2

% u_x
D1u = zeros(size(u));
for n = -p:p
    D1u = D1u + d1(n+p+1).*circshift(u,-n);
end
% c*u_x
D1u = c.*D1u;
% (c*u_x)_x
res = zeros(size(u));
for n = -p:p
    res = res + d1(n+p+1).*circshift(D1u,-n);
end
% b*(c*u_x)_x + F
res = b.*res + F(t);

end
