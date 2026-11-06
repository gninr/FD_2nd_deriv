% narrow-stencil finite difference weights

%% setup
% order of accuracy 2p
p = 5;
disp("Order")
disp(2*p)

%% assemble matrix
M = sym(zeros(2*p+1));
for k = 0:p
    d = nchoosek(p,k)^2/nchoosek(2*p,p);
    a = compute_a(p,k);
    M = M + d*(a*a');
end
disp("M =")
disp(M)

function coeff = compute_a(m,k)
coeff = sym(zeros(2*m+1,1));
n = sym(k-m:k);
coeff(n+m+1) = (-1).^(n+1).*factorial(k).*factorial(m-k)./(n.*factorial(k-n).*factorial(m-k+n));
coeff(m+1) = harmonic(sym(m-k)) - harmonic(sym(k));
end

%% compute weights
B = sym(zeros(2*p+1));
for m = -p:p
    for n = -p:p
        if abs(m-n) > p
            continue
        end
        B(m+p+1,n+p+1) = -M(n-m+p+1,-m+p+1);
    end
end
disp("B =")
disp(B)

%% output
disp(latex(B))
