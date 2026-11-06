function res = dvdt_narrow(u,t,b,c,F,B,p)
% Approximate v_t = b*(c*u_x)_x+F using 2p-th order narrow stencil B

% (c*u_x)_x
res = zeros(size(u));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        res = res + B(m+p+1,n+p+1)*circshift(c,-m).*circshift(u,-n);
    end
end
% b*(c*u_x)_x + F
res = b.*res + F(t);

end
