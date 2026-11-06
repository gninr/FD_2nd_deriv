% 1D pi-mode test
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
T = 1;
% coordinate transformation
epsilon = 0.1;
x = @(xi) xi+epsilon*sin(2*pi*xi);
dx_dxi = @(xi) 1+epsilon*2*pi*cos(2*pi*xi);
% initial conditions
f = @(x) exp(-500*(x-0.5).^2);
g = @(x) zeros(size(x));
% exact solution
sol = @(x,t) (f(mod(x+t,1))+f(mod(x-t,1)))/2;

%% discretization
% temporal discretization
dt = 1e-4;

% spatial discretization
N = 64;
h = 1/N;
xi = (1:N)*h;
c = 1./dx_dxi(xi);
Fh = @(t) zeros(size(xi));
u0 = f(x(xi));
v0 = g(x(xi));
ref = sol(x(xi),T);

%% narrow stencil
B = double(B1/h^2);
u_narrow = rk4(@(u,t)dvdt_narrow(u,t,c,c,Fh,B,p),u0,v0,T,dt);

%% wide stencil
d1 = double(D/h);
u_wide = rk4(@(u,t)dvdt_wide(u,t,c,c,Fh,d1,p),u0,v0,T,dt);

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    u_kp = rk4(@(u,t)dvdt_narrow(u,t,c,c,Fh,B,p),u0,v0,T,dt);
end

%% plot error
if p == 3
    subplot(2,1,1)
    plot(xi,u_narrow,xi,u_wide,xi,u_kp,xi,ref,'k--')
    xlabel("\xi")
    ylim([-1 1])
    legend("Narrow","Wide","KP","Exact")
    title("T = "+num2str(T))

    subplot(2,1,2)
    plot(xi,u_narrow-ref,xi,u_wide-ref,xi,u_kp-ref)
    xlabel("\xi")
    legend("Narrow","Wide","KP")
    title("Difference")
else
    subplot(2,1,1)
    plot(xi,u_narrow,xi,u_wide,xi,ref,'k--')
    xlabel("\xi")
    ylim([-1 1])
    legend("Narrow","Wide","Exact")
    title("T = "+num2str(T))

    subplot(2,1,2)
    plot(xi,u_narrow-ref,xi,u_wide-ref)
    xlabel("\xi")
    legend("Narrow","Wide")
    title("Difference")
end
sgtitle("Order "+int2str(2*p))

%% output
folder = "../../data/scalar1/";
% difference
print_u(xi,u_narrow,folder,"narrow"+int2str(p));
print_u(xi,u_wide,folder,"wide"+int2str(p));
if p == 3
    print_u(xi,u_kp,folder,"kp");
end

function print_u(xi,u,folder,fn)
    fname = folder+"pi_"+fn+".dat";
    fid = fopen(fname,'w+');
    fprintf(fid,"xi u\n"+sprintf("%f %f\n",[0 xi;u(end) u]));
    fclose(fid);
end
