% 2D accuracy test
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
x = @(x1,x2) x1+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
y = @(x1,x2) x2+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
dx_dx1 = @(x1,x2) 1+epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dx_dx2 = @(x1,x2)   epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
dy_dx1 = @(x1,x2)   epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dy_dx2 = @(x1,x2) 1+epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
% exact solution
sol = @(x,y,t) sin(2*pi*x).*sin(2*pi*y)*cos(2*sqrt(2)*pi*t);
% initial condition
f = @(x,y) sin(2*pi*x).*sin(2*pi*y);
g = @(x,y) zeros(size(x));

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
    x_ = (1:N)*h;
    [x1,x2] = meshgrid(x_,x_);
    J11 = dx_dx1(x1,x2);
    J12 = dx_dx2(x1,x2);
    J21 = dy_dx1(x1,x2);
    J22 = dy_dx2(x1,x2);
    J = J11.*J22-J12.*J21;
    b = 1./J;
    c11 = b.*(J12.^2+J22.^2);
    c12 = -b.*(J11.*J12+J21.*J22);
    c22 = b.*(J11.^2+J21.^2);
    Fh = @(t) zeros(size(x1)); % forcing
    u0 = f(x(x1,x2),y(x1,x2));
    v0 = g(x(x1,x2),y(x1,x2));
    ref = sol(x(x1,x2),y(x1,x2),T);

    d1 = double(D/h);

    %%% narrow stencil
    B = double(B1/h^2);
    u = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T,dt);
    err_narrow(m) = sqrt(h^2*sum((u-ref).^2,'all'));

    %%% wide stencil
    u = rk4(@(u,t)dvdt_wide(u,t,b,c11,c12,c22,Fh,d1,p),u0,v0,T,dt);
    err_wide(m) = sqrt(h^2*sum((u-ref).^2,'all'));

    %%% Kamakoti & Pantano
    if p == 3
        B = double(B2/h^2);
        u = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T,dt);
        err_kp(m) = sqrt(h^2*sum((u-ref).^2,'all'));
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
s = sprintf("$%d^2$",Ns(1));
s = add(s,err_narrow(1));
s = add(s,err_wide(1));
if p == 3
    s = add(s,err_kp(1));
end
s = s + " \\\\\n";
% rest of rows
add = @(s,err,rate) s + sprintf(" & %.2e & %.2f",err,rate);
for m = 2:ntest
    s = s + sprintf("$%d^2$",Ns(m));
    s = add(s,err_narrow(m),rate_narrow(m-1));
    s = add(s,err_wide(m),rate_wide(m-1));
    if p == 3
        s = add(s,err_kp(m),rate_kp(m-1));
    end
    s = s + " \\\\\n";
end
fprintf(s)
