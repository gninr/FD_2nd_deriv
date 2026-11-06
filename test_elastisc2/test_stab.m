% elastic stability test
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
T = 50;
T_kp = min(T,5);
% coordinate transformation
epsilon1 = 0.1;
x = @(x1,x2) x1+epsilon1*sin(2*pi*x1).*sin(2*pi*x2);
y = @(x1,x2) x2+epsilon1*sin(2*pi*x1).*sin(2*pi*x2);
dx_dx1 = @(x1,x2) 1+epsilon1*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dx_dx2 = @(x1,x2)   epsilon1*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
dy_dx1 = @(x1,x2)   epsilon1*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dy_dx2 = @(x1,x2) 1+epsilon1*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
% coefficients
epsilon2 = 0.8;
k = 20*pi;
rho    = @(x,y) 1+epsilon2*cos(k*x).*sin(k*y);
lambda = @(x,y) 1+epsilon2*sin(k*x).*sin(k*y);
mu     = @(x,y) 1+epsilon2*cos(k*x).*cos(k*y);
% exact solution
sol  = @(x,y,t) cat(3, ...
    sin(2*pi*(x+t)).*sin(2*pi*y), ...
    sin(2*pi*(x+2*t)).*sin(2*pi*y));
% auxiliary functions
ut   = @(x,y,t) cat(3, ...
    2*pi*cos(2*pi*(x+t)).*sin(2*pi*y), ...
    4*pi*cos(2*pi*(x+2*t)).*sin(2*pi*y));
utt  = @(x,y,t) cat(3, ...
    -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y), ...
    -16*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y));
u1x  = @(x,y,t) 2*pi*cos(2*pi*(x+t)).*sin(2*pi*y);
u2x  = @(x,y,t) 2*pi*cos(2*pi*(x+2*t)).*sin(2*pi*y);
u1y  = @(x,y,t) 2*pi*sin(2*pi*(x+t)).*cos(2*pi*y);
u2y  = @(x,y,t) 2*pi*sin(2*pi*(x+2*t)).*cos(2*pi*y);
u1xx = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
u2xx = @(x,y,t) -4*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y);
u1xy = @(x,y,t) 4*pi^2*cos(2*pi*(x+t)).*cos(2*pi*y);
u2xy = @(x,y,t) 4*pi^2*cos(2*pi*(x+2*t)).*cos(2*pi*y);
u1yy = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
u2yy = @(x,y,t) -4*pi^2*sin(2*pi*(x+2*t)).*sin(2*pi*y);
lambdax = @(x,y) k*epsilon2*cos(k*x).*sin(k*y);
lambday = @(x,y) k*epsilon2*sin(k*x).*cos(k*y);
mux =     @(x,y) -k*epsilon2*sin(k*x).*cos(k*y);
muy =     @(x,y) -k*epsilon2*cos(k*x).*sin(k*y);
% forcing
F = @(x,y,t) rho(x,y).*utt(x,y,t)-cat(3, ...
    (lambdax(x,y)+2*mux(x,y)).*u1x(x,y,t)+lambdax(x,y).*u2y(x,y,t)+(lambda(x,y)+2*mu(x,y)).*u1xx(x,y,t)+lambda(x,y).*u2xy(x,y,t) ...
    +muy(x,y).*(u2x(x,y,t)+u1y(x,y,t))+mu(x,y).*(u2xy(x,y,t)+u1yy(x,y,t)), ...
    mux(x,y).*(u2x(x,y,t)+u1y(x,y,t))+mu(x,y).*(u2xx(x,y,t)+u1xy(x,y,t)) ...
    +lambday(x,y).*u1x(x,y,t)+(lambday(x,y)+2*muy(x,y)).*u2y(x,y,t)+lambda(x,y).*u1xy(x,y,t)+(lambda(x,y)+2*mu(x,y)).*u2yy(x,y,t));
% initial conditions
f = @(x,y) sol(x,y,0);
g = @(x,y) ut(x,y,0);

%% discretization
% temporal discretization
dt = 1e-4;

% spatial discretization
N = 64;
h = 1/N;
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
b = 1./rho(x(x1,x2),y(x1,x2))./J;
muh = mu(x(x1,x2),y(x1,x2));
lambdah = lambda(x(x1,x2),y(x1,x2));
kappa1 = @(k,l) J.*((lambdah+2*muh).*Jinv{k,1}.*Jinv{l,1}+muh.*Jinv{k,2}.*Jinv{l,2});
kappa2 = @(k,l) J.*(lambdah.*Jinv{k,1}.*Jinv{l,2}+muh.*Jinv{k,2}.*Jinv{l,1});
kappa3 = @(k,l) J.*(muh.*Jinv{k,1}.*Jinv{l,1}+(lambdah+2*muh).*Jinv{k,2}.*Jinv{l,2});
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
Fh = @(t) F(x(x1,x2),y(x1,x2),t)./rho(x(x1,x2),y(x1,x2));
u0 = f(x(x1,x2),y(x1,x2));
v0 = g(x(x1,x2),y(x1,x2));
ref = @(t) sol(x(x1,x2),y(x1,x2),t);
d1 = double(D/h);

%% narrow stencil
B = double(B1/h^2);
[~,t,err_narrow] = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T,dt,ref);
err_narrow = h*err_narrow;

%% wide stencil
[~,~,err_wide] = rk4(@(u,t)dvdt_wide(u,t,b,c,Fh,d1,p),u0,v0,T,dt,ref);
err_wide = h*err_wide;

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    [~,t_kp,err_kp] = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T_kp,dt,ref);
    err_kp = h*err_kp;
end

%% plot
if p == 3
    semilogy(t,err_narrow,t,err_wide,t_kp,err_kp)
    xlabel("t")
    legend("Narrow","Wide","KP")
else
    semilogy(t,err_narrow,t,err_wide)
    xlabel("t")
    legend("Narrow","Wide")
end
title("Order "+int2str(2*p))

%% output
folder = "../../data/elastic2/";
% error
print_err(t,err_narrow,folder,"narrow"+int2str(p));
print_err(t,err_wide,folder,"wide"+int2str(p));
if p == 3
    print_err(t_kp,err_kp,folder,"kp");
end

function print_err(t,err,folder,fn)
    fname = folder+"evol_"+fn+".dat";
    fid = fopen(fname,'w+');
    fprintf(fid,"t err\n"+sprintf("%f %f\n",[t(1:100:end);err(1:100:end)]));
    fclose(fid);
end
