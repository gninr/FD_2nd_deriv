% 1D stability test
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
rho = @(x) 1+epsilon*cos(k*x);
kappa = @(x) 1+epsilon*sin(k*x);
% exact solution
sol = @(x,t) sin(2*pi*(x+t));
% auxiliary functions
ut     = @(x,t) 2*pi*cos(2*pi*(x+t));
utt    = @(x,t) -4*pi^2*sin(2*pi*(x+t));
ux     = @(x,t) 2*pi*cos(2*pi*(x+t));
uxx    = @(x,t) -4*pi^2*sin(2*pi*(x+t));
kappax = @(x)   k*epsilon*cos(k*x);
% forcing
F = @(x,t) rho(x).*utt(x,t)-kappax(x).*ux(x,t)-kappa(x).*uxx(x,t);
% initial conditions
f = @(x) sol(x,0);
g = @(x) ut(x,0);

%% discretization
% temporal discretization
dt = 1e-4;

% spatial discretization
N = 64;
h = 1/N;
x = (1:N)*h;
b = 1./rho(x);
c = kappa(x);
Fh = @(t) F(x,t)./rho(x);
u0 = f(x);
v0 = g(x);
ref = @(t) sol(x,t);

%% narrow stencil
B = double(B1/h^2);
[~,t,err_narrow] = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,B,p),u0,v0,T,dt,ref);
err_narrow = sqrt(h)*err_narrow;

M = zeros(2*p+1);
for m = 1:2*p+1
    M(m,:) = -circshift(B(end-m+1,:),m-p-1);
end
D2_elem = zeros(N);
D2_elem(1:2*p+1,1:2*p+1) = -M;
D2_elem = circshift(D2_elem,[-p -p]);
D2 = zeros(N);
for j = 1:N
    D2 = D2 + c(j)*circshift(D2_elem,[j-1 j-1]);
end
D2 = diag(b)*D2;
eig_narrow = eig(D2);

%% wide stencil
d1 = double(D/h);
[~,~,err_wide] = rk4(@(u,t)dvdt_wide(u,t,b,c,Fh,d1,p),u0,v0,T,dt,ref);
err_wide = sqrt(h)*err_wide;

D1_elem = zeros(1,N);
D1_elem(1:2*p+1) = d1;
D1_elem = circshift(D1_elem,-p);
D1 = zeros(N);
for j = 1:N
    D1(j,:) = circshift(D1_elem,j-1);
end
D2 = diag(b)*D1*diag(c)*D1;
eig_wide = eig(D2);

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    [~,t_kp,err_kp] = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,B,p),u0,v0,T_kp,dt,ref);
    err_kp = sqrt(h)*err_kp;

    M = zeros(2*p+1);
    for m = 1:2*p+1
        M(m,:) = -circshift(B(end-m+1,:),m-p-1);
    end
    D2_elem = zeros(N);
    D2_elem(1:2*p+1,1:2*p+1) = -M;
    D2_elem = circshift(D2_elem,[-p -p]);
    D2 = zeros(N);
    for j = 1:N
        D2 = D2 + c(j)*circshift(D2_elem,[j-1 j-1]);
    end
    D2 = diag(b)*D2;
    eig_kp = eig(D2);
end

%% plot error
if p == 3
    subplot(2,1,1)
    semilogy(t,err_narrow,t,err_wide,t_kp,err_kp)
    xlabel("t")
    legend("Narrow","Wide","KP")

    subplot(2,3,4)
    scatter(real(eig_narrow),imag(eig_narrow))
    title("Narrow")

    subplot(2,3,5)
    scatter(real(eig_wide),imag(eig_wide))
    title("Wide")

    subplot(2,3,6)
    scatter(real(eig_kp),imag(eig_kp))
    title("KP")
else
    subplot(2,1,1)
    semilogy(t,err_narrow,t,err_wide)
    xlabel("t")
    legend("Narrow","Wide")

    subplot(2,2,3)
    scatter(real(eig_narrow),imag(eig_narrow))
    title("Narrow")

    subplot(2,2,4)
    scatter(real(eig_wide),imag(eig_wide))
    title("Wide")
end
sgtitle("Order "+int2str(2*p))

%% output
folder = "../../data/scalar1/";
% error
print_err(t,err_narrow,folder,"narrow"+int2str(p));
print_err(t,err_wide,folder,"wide"+int2str(p));
if p == 3
    print_err(t_kp,err_kp,folder,"kp");
end
% eigenvalues
print_eig(eig_narrow,folder,"narrow"+int2str(p));
print_eig(eig_wide,folder,"wide"+int2str(p));
if p == 3
    print_eig(eig_kp,folder,"kp");
end

function print_err(t,err,folder,fn)
    fname = folder+"evol_"+fn+".dat";
    fid = fopen(fname,'w+');
    fprintf(fid,"t err\n"+sprintf("%f %f\n",[t(1:100:end);err(1:100:end)]));
    fclose(fid);
end

function print_eig(eigval,folder,fn)
    fname = folder+"eig_"+fn+".dat";
    fid = fopen(fname,'w+');
    fprintf(fid,"real imag\n"+sprintf("%f %f\n",[real(eigval) imag(eigval)]'));
    fclose(fid);
end
