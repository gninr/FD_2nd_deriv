% 2D stability test
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
% coefficients
epsilon = 0.8;
k = 20*pi;
rho = @(x,y) 1+epsilon*cos(k*x).*sin(k*y);
kappa11 = @(x,y) 1+epsilon*sin(k*x).*sin(k*y);
kappa12 = @(x,y)   epsilon*sin(k*x).*cos(k*y);
kappa22 = @(x,y) 1+epsilon*cos(k*x).*cos(k*y);
% exact solution
sol = @(x,y,t) sin(2*pi*(x+t)).*sin(2*pi*y);
% auxiliary functions
ut       = @(x,y,t) 2*pi*cos(2*pi*(x+t)).*sin(2*pi*y);
utt      = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
ux       = @(x,y,t) 2*pi*cos(2*pi*(x+t)).*sin(2*pi*y);
uy       = @(x,y,t) 2*pi*sin(2*pi*(x+t)).*cos(2*pi*y);
uxx      = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
uxy      = @(x,y,t) 4*pi^2*cos(2*pi*(x+t)).*cos(2*pi*y);
uyy      = @(x,y,t) -4*pi^2*sin(2*pi*(x+t)).*sin(2*pi*y);
kappa11x = @(x,y)   k*epsilon*cos(k*x).*sin(k*y);
kappa12x = @(x,y)   k*epsilon*cos(k*x).*cos(k*y);
kappa12y = @(x,y)   -k*epsilon*sin(k*x).*sin(k*y);
kappa22y = @(x,y)   -k*epsilon*cos(k*x).*sin(k*y);
% forcing
F = @(x,y,t) rho(x,y).*utt(x,y,t) ...
    -kappa11x(x,y).*ux(x,y,t)-kappa11(x,y).*uxx(x,y,t) ...
    -kappa12x(x,y).*uy(x,y,t)-2*kappa12(x,y).*uxy(x,y,t)-kappa12y(x,y).*ux(x,y,t) ...
    -kappa22y(x,y).*uy(x,y,t)-kappa22(x,y).*uyy(x,y,t);
% initial condition
f = @(x,y) sol(x,y,0);
g = @(x,y) ut(x,y,0);

%% discretization
% temporal discretization
dt = 1e-4;

% spatial discretization
N = 64;
h = 1/N;
x_ = (1:N)*h;
[x,y] = meshgrid(x_,x_);
b = 1./rho(x,y);
c11 = kappa11(x,y);
c12 = kappa12(x,y);
c22 = kappa22(x,y);
Fh = @(t) F(x,y,t)./rho(x,y);
u0 = f(x,y);
v0 = g(x,y);
ref = @(t) sol(x,y,t);
d1 = double(D/h);

%% narrow stencil
B = double(B1/h^2);
[~,t,err_narrow] = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T,dt,ref);
err_narrow = h*err_narrow;

%% wide stencil
[~,~,err_wide] = rk4(@(u,t)dvdt_wide(u,t,b,c11,c12,c22,Fh,d1,p),u0,v0,T,dt,ref);
err_wide = h*err_wide;

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    [~,t_kp,err_kp] = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T_kp,dt,ref);
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
folder = "../../data/scalar2/";
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
