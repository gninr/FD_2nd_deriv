function B = int_narrow(p)
% 2p-th order narrow-stencil operator

% compute M
M = sym(zeros(2*p+1));
for n = 0:p
    d = nchoosek(p,n)^2/nchoosek(2*p,p);
    a = compute_a(p,n);
    M = M + d*(a*a');
end

% compute B
B = sym(zeros(2*p+1));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        B(m+p+1,n+p+1) = -M(n-m+p+1,-m+p+1);
    end
end
end

function coeff = compute_a(s,n)
coeff = sym(zeros(2*s+1,1));
k = sym(n-s:n);
coeff(k+s+1) = (-1).^(k+1).*factorial(n).*factorial(s-n)./(k.*factorial(n-k).*factorial(s-n+k));
coeff(s+1) = harmonic(sym(s-n)) - harmonic(sym(n));
end
