% 2D pi-mode test
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
x = @(x1,x2) x1+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
y = @(x1,x2) x2+epsilon*sin(2*pi*x1).*sin(2*pi*x2);
dx_dx1 = @(x1,x2) 1+epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dx_dx2 = @(x1,x2)   epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
dy_dx1 = @(x1,x2)   epsilon*2*pi*cos(2*pi*x1).*sin(2*pi*x2);
dy_dx2 = @(x1,x2) 1+epsilon*2*pi*sin(2*pi*x1).*cos(2*pi*x2);
% exact solution
sx = 250;
sy = 500;
sol = @(x,y,t) exp(-sx*(mod(x+t,1)-0.5).^2-sy*(mod(y+t,1)-0.5).^2);
% auxiliary functions
ut  = @(x,y,t) (-2*sx*(mod(x+t,1)-0.5)-2*sy*(mod(y+t,1)-0.5)).*sol(x,y,t);
utt = @(x,y,t) (-2*sx-2*sy+(-2*sx*(mod(x+t,1)-0.5)-2*sy*(mod(y+t,1)-0.5)).^2).*sol(x,y,t);
uxx = @(x,y,t) (-2*sx+(-2*sx*(mod(x+t,1)-0.5)).^2).*sol(x,y,t);
uyy = @(x,y,t) (-2*sy+(-2*sy*(mod(y+t,1)-0.5)).^2).*sol(x,y,t);
% forcing
F = @(x,y,t) utt(x,y,t)-uxx(x,y,t)-uyy(x,y,t);
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
b = 1./J;
c11 = b.*(J12.^2+J22.^2);
c12 = -b.*(J11.*J12+J21.*J22);
c22 = b.*(J11.^2+J21.^2);
Fh = @(t) F(x(x1,x2),y(x1,x2),t); % forcing
u0 = f(x(x1,x2),y(x1,x2));
v0 = g(x(x1,x2),y(x1,x2));

d1 = double(D/h);

%% narrow stencil
B = double(B1/h^2);
u_narrow = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T,dt);

%% wide stencil
u_wide = rk4(@(u,t)dvdt_wide(u,t,b,c11,c12,c22,Fh,d1,p),u0,v0,T,dt);

%% Kamakoti & Pantano
if p == 3
    B = double(B2/h^2);
    u_kp = rk4(@(u,t)dvdt_narrow(u,t,b,c11,c12,c22,Fh,d1,B,p),u0,v0,T,dt);
end

%% plot
x_ref = (1:4*N)*h/4;
[x1_ref,x2_ref] = meshgrid(x_ref,x_ref);
ref = sol(x(x1_ref,x2_ref),y(x1_ref,x2_ref),T);
cmap = flipud(colormap('gray'));
colormap(cmap)
cl = [-0.05,0.05];
if p == 3
    subplot(2,2,1)
    imagesc(x_,x_,u_narrow)
    colorbar
    clim(cl)
    title("Narrow")

    subplot(2,2,2)
    imagesc(x_,x_,u_wide)
    colorbar
    clim(cl)
    title("Wide")

    subplot(2,2,3)
    imagesc(x_,x_,u_kp)
    colorbar
    clim(cl)
    title("KP")

    subplot(2,2,4)
    imagesc(x_ref,x_ref,ref)
    colorbar
    clim(cl)
    title("Exact")
else
    subplot(2,2,1)
    imagesc(x_,x_,u_narrow)
    colorbar
    clim(cl)
    title("Narrow")

    subplot(2,2,2)
    imagesc(x_,x_,u_wide)
    colorbar
    clim(cl)
    title("Wide")

    subplot(2,2,4)
    imagesc(x_ref,x_ref,ref)
    colorbar
    clim(cl)
    title("Exact")
end
sgtitle("Order "+int2str(2*p)+", T = "+num2str(T))

%% output
folder = "../../data/scalar2/";
% colormap
str_map = "\\pgfplotsset{\n\tcolormap={matlab}{\n"+sprintf("\t\trgb=(%g,%g,%g)\n",cmap')+sprintf("\t},\n\tpoint meta min=%f,\n\tpoint meta max=%f,\n}\n",cl);
fid = fopen(folder+"pi_map.tex",'w+');
fprintf(fid,str_map);
fclose(fid);
% solution
ncolor = size(cmap,1);
conv_img = @(u) floor((flip(u)-cl(1))/(cl(2)-cl(1))*(ncolor-1))+1;
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
