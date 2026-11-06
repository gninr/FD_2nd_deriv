% elastic accuracy test
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
sol = @(x,y,t) cat(3, ...
    sin(2*pi*(x+t)).*sin(2*pi*y), ...
    sin(2*pi*(x+2*t)).*sin(2*pi*y));
% auxiliary functions
ut  = @(x,y,t) cat(3, ...
    2*pi*cos(2*pi*(x+t)).*sin(2*pi*y), ...
    4*pi*cos(2*pi*(x+2*t)).*sin(2*pi*y));
utt = @(x,y,t) cat(3, ...
    -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y), ...
    -16*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y));
u1xx = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
u2xx = @(x,y,t) -4*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y);
u1xy = @(x,y,t) 4*pi^2*cos(2*pi*(x+t)).*cos(2*pi*y);
u2xy = @(x,y,t) 4*pi^2*cos(2*pi*(x+2*t)).*cos(2*pi*y);
u1yy = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
u2yy = @(x,y,t) -4*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y);
% forcing
F = @(x,y,t) utt(x,y,t)-cat(3, ...
    3*u1xx(x,y,t)+u1yy(x,y,t)+2*u2xy(x,y,t), ...
    2*u1xy(x,y,t)+u2xx(x,y,t)+3*u2yy(x,y,t));
% initial conditions
f = @(x,y) sol(x,y,0);
g = @(x,y) ut(x,y,0);

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
    Jinv{1,1} = J22./J;
    Jinv{1,2} = -J12./J;
    Jinv{2,1} = -J21./J;
    Jinv{2,2} = J11./J;
    b = 1./J;
    kappa1 = @(k,l) J.*(3*Jinv{k,1}.*Jinv{l,1}+Jinv{k,2}.*Jinv{l,2});
    kappa2 = @(k,l) J.*(Jinv{k,1}.*Jinv{l,2}+Jinv{k,2}.*Jinv{l,1});
    kappa3 = @(k,l) J.*(Jinv{k,1}.*Jinv{l,1}+3*Jinv{k,2}.*Jinv{l,2});
    c = cell(2,2,2); % d(c(k,l,i)*du_i/dx_l)/dx_k
    % u1
    c{1,1,1} = cat(3,kappa1(1,1),kappa2(1,1));
    c{1,2,1} = cat(3,kappa1(1,2),kappa2(2,1));
    c{2,1,1} = cat(3,kappa1(2,1),kappa2(1,2));
    c{2,2,1} = cat(3,kappa1(2,2),kappa2(2,2));
    % u2
    c{1,1,2} = cat(3,kappa2(1,1),kappa3(1,1));
    c{1,2,2} = cat(3,kappa2(1,2),kappa3(1,2));
    c{2,1,2} = cat(3,kappa2(2,1),kappa3(2,1));
    c{2,2,2} = cat(3,kappa2(2,2),kappa3(2,2));
    Fh = @(t) F(x(x1,x2),y(x1,x2),t);
    u0 = f(x(x1,x2),y(x1,x2));
    v0 = g(x(x1,x2),y(x1,x2));
    ref = sol(x(x1,x2),y(x1,x2),T);

    d1 = double(D/h);

    %%% narrow stencil
    B = double(B1/h^2);
    u = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T,dt);
    err_narrow(m) = sqrt(h^2*sum((u-ref).^2,'all'));

    %%% wide stencil
    u = rk4(@(u,t)dvdt_wide(u,t,b,c,Fh,d1,p),u0,v0,T,dt);
    err_wide(m) = sqrt(h^2*sum((u-ref).^2,'all'));

    %%% Kamakoti & Pantano
    if p == 3
        B = double(B2/h^2);
        u = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T,dt);
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
