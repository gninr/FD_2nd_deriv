% elastic pi-mode test
addpath("../util/")

%% stencil
% order of accuracy 2p
p = 5;
D = int_d1(p);
B1 = int_narrow(p);
if p == 3
    B2 = int_kp;
end

%% problem setup
T = 1;
% coordinate transformation
epsilon = 0.1;
x = @(x1,x2) x1+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
y = @(x1,x2) x2+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
dx_dx1 = @(x1,x2) 1+epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dx_dx2 = @(x1,x2)   epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
dy_dx1 = @(x1,x2)   epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dy_dx2 = @(x1,x2) 1+epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
% exact solution
sx = 250;
sy = 500;
sol1 = @(x,y,t) exp(-sx*(mod(x+t,1)-0.5).^2-sy*(mod(y+t,1)-0.5).^2);
sol2 = @(x,y,t) exp(-sx*(mod(x+2*t,1)-0.5).^2-sy*(mod(y+t,1)-0.5).^2);
sol = @(x,y,t) cat(3,sol1(x,y,t),sol2(x,y,t));
% auxiliary functions
ut  = @(x,y,t) sol(x,y,t).*cat(3, ...
    -2*sx*(mod(x+t,1)-0.5)-2*sy*(mod(y+t,1)-0.5), ...
    -4*sx*(mod(x+2*t,1)-0.5)-2*sy*(mod(y+t,1)-0.5));
utt = @(x,y,t) sol(x,y,t).*cat(3, ...
    -2*sx-2*sy+(-2*sx*(mod(x+t,1)-0.5)-2*sy*(mod(y+t,1)-0.5)).^2, ...
    -8*sx-2*sy+(-4*sx*(mod(x+2*t,1)-0.5)-2*sy*(mod(y+t,1)-0.5)).^2);
u1xx = @(x,y,t) (-2*sx+(-2*sx*(mod(x+t,1)-0.5)).^2)                 .*sol1(x,y,t);
u2xx = @(x,y,t) (-2*sx+(-2*sx*(mod(x+2*t,1)-0.5)).^2)               .*sol2(x,y,t);
u1xy = @(x,y,t) (-2*sx*(mod(x+t,1)-0.5)).*(-2*sy*(mod(y+t,1)-0.5))  .*sol1(x,y,t);
u2xy = @(x,y,t) (-2*sx*(mod(x+2*t,1)-0.5)).*(-2*sy*(mod(y+t,1)-0.5)).*sol2(x,y,t);
u1yy = @(x,y,t) (-2*sy+(-2*sy*(mod(y+t,1)-0.5)).^2)                 .*sol1(x,y,t);
u2yy = @(x,y,t) (-2*sy+(-2*sy*(mod(y+t,1)-0.5)).^2)                 .*sol2(x,y,t);
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

d1 = double(D/h);

%% narrow stencil
B = double(B1/h^2);
u_narrow = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T,dt);

%% wide stencil
u_wide = rk4(@(u,t)dvdt_wide(u,t,b,c,Fh,d1,p),u0,v0,T,dt);

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    u_kp = rk4(@(u,t)dvdt_narrow(u,t,b,c,Fh,d1,B,p),u0,v0,T,dt);
end

%% plot
magn = @(u) sqrt(u(:,:,1).^2+u(:,:,2).^2);
x_ref = (1:4*N)*h/4;
[x1_ref,x2_ref] = meshgrid(x_ref,x_ref);
ref = sol(x(x1_ref,x2_ref),y(x1_ref,x2_ref),T);
cmap = flipud(colormap('gray'));
colormap(cmap)
cl = [0,0.09];
if p == 3
    subplot(2,2,1)
    imagesc(magn(u_narrow))
    colorbar
    clim(cl)
    title("Narrow")

    subplot(2,2,2)
    imagesc(magn(u_wide))
    colorbar
    clim(cl)
    title("Wide")

    subplot(2,2,3)
    imagesc(magn(u_kp))
    colorbar
    clim(cl)
    title("KP")

    subplot(2,2,4)
    imagesc(magn(ref))
    colorbar
    clim(cl)
    title("Exact")
else
    subplot(2,2,1)
    imagesc(magn(u_narrow))
    colorbar
    clim(cl)
    title("Narrow")

    subplot(2,2,2)
    imagesc(magn(u_wide))
    colorbar
    clim(cl)
    title("Wide")

    subplot(2,2,4)
    imagesc(magn(ref))
    colorbar
    clim(cl)
    title("Exact")
end
sgtitle("Order "+int2str(2*p)+", T = "+num2str(T))

%% output
folder = "../../data/elastic2/";
% colormap
str_map = "\\pgfplotsset{\n\tcolormap={matlab}{\n"+sprintf("\t\trgb=(%g,%g,%g)\n",cmap')+sprintf("\t},\n\tpoint meta min=%f,\n\tpoint meta max=%f,\n}\n",cl);
fid = fopen(folder+"pi_map.tex",'w+');
fprintf(fid,str_map);
fclose(fid);
% solution
prefix = folder+"pi_";
ncolor = size(cmap,1);
conv_img = @(u) floor((magn(flip(u))-cl(1))/(cl(2)-cl(1))*(ncolor-1))+1;
figure
print_u(conv_img(ref),cmap,folder,"ref");
print_u(conv_img(u_narrow),cmap,folder,"narrow"+int2str(p));
print_u(conv_img(u_wide),cmap,folder,"wide"+int2str(p));
if p == 3
    print_u(conv_img(u_kp),cmap,folder,"kp");
end

function print_u(img,cmap,folder,fn)
    fname = folder+"pi_"+fn+".eps";
    imshow(img,cmap)
    saveas(gca,fname)
end
