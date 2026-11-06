% 1D accuracy test
addpath("../util/")

%% stencil
% order of accuracy 2p
p = 3;
D = int_d1(p);
B1 = int_narrow(p);
if p == 3
    B2 = int_kp;
end

%% problem setup
T = 0.5;
% coordinate transformation
epsilon = 0.1;
x = @(xi) xi+epsilon*sin(2*pi*xi);
dx_dxi = @(xi) 1+epsilon*2*pi*cos(2*pi*xi);
% exact solution
sol = @(x,t) sin(2*pi*x)*cos(2*pi*t);
% initial conditions
f = @(x) sin(2*pi*x);
g = @(x) zeros(size(x));

%% discretization
% temporal discretization
dt = 1e-4;

% spatial discretization
Ns = 16*2.^(0:3);
hs = 1./Ns;

% error
ntest = length(Ns);
err_narrow = zeros(1,ntest);
err_wide = zeros(1,ntest);
err_kp = zeros(1,ntest);

for m = 1:ntest
    N = Ns(m);
    h = hs(m);
    xi = (1:N)*h;
    c = 1./dx_dxi(xi);
    Fh = @(t) zeros(size(x)); % forcing
    u0 = f(x(xi));
    v0 = g(x(xi));
    ref = sol(x(xi),T);

    %%% narrow stencil
    B = double(B1/h^2);
    u = rk4(@(u,t)dvdt_narrow(u,t,c,c,Fh,B,p),u0,v0,T,dt);
    err_narrow(m) = sqrt(h*sum((u-ref).^2));

    %%% wide stencil
    d1 = double(D/h);
    u = rk4(@(u,t)dvdt_wide(u,t,c,c,Fh,d1,p),u0,v0,T,dt);
    err_wide(m) = sqrt(h*sum((u-ref).^2));

    %%% Kamakoti & Pantano
    if p == 3
        B = double(B2/h^2);
        u = rk4(@(u,t)dvdt_narrow(u,t,c,c,Fh,B,p),u0,v0,T,dt);
        err_kp(m) = sqrt(h*sum((u-ref).^2));
    end
end

%% plot error
hs_2p = hs.^(2*p);
loglog(hs,err_narrow,hs,err_wide,hs,err_kp,hs,hs_2p,'k--')
xlabel("h")
legend("Narrow","Wide","KP","O(h^{2p})")
title("Order "+int2str(2*p)+", l_2 error")

%% output
% convergence rate
rate_narrow = -diff(log(err_narrow))/log(2);
rate_wide = -diff(log(err_wide))/log(2);
rate_kp = -diff(log(err_kp))/log(2);
% first row
add = @(s,err) s + sprintf(" & %.2e & ",err);
s = sprintf("%d",Ns(1));
s = add(s,err_narrow(1));
s = add(s,err_wide(1));
if p == 3
    s = add(s,err_kp(1));
end
s = s + " \\\\\n";
% rest of rows
add = @(s,err,rate) s + sprintf(" & %.2e & %.2f",err,rate);
for m = 2:ntest
    s = s + sprintf("%d",Ns(m));
    s = add(s,err_narrow(m),rate_narrow(m-1));
    s = add(s,err_wide(m),rate_wide(m-1));
    if p == 3
        s = add(s,err_kp(m),rate_kp(m-1));
    end
    s = s + " \\\\\n";
end
fprintf(s)
