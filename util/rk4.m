function [u,t,err] = rk4(dvdt,u,v,T,dt,ref)
% Fourth order Runge Kutta for u'=v, v'=dvdt(u)

t = 0:dt:T;
Nt = length(t);
if nargin == 6
    err = zeros(1,Nt);
end
for n = 1:Nt-1
    tn = t(n);
    if nargin == 6
        err(n) = norm(ref(tn)-u,'fro');
    end

    dudt1 = v;
    dvdt1 = dvdt(u,tn);

    dudt2 = v+dt/2*dvdt1;
    dvdt2 = dvdt(u+dt/2*dudt1,tn+dt/2);

    dudt3 = v+dt/2*dvdt2;
    dvdt3 = dvdt(u+dt/2*dudt2,tn+dt/2);

    dudt4 = v+dt*dvdt3;
    dvdt4 = dvdt(u+dt*dudt3,tn+dt);

    u = u+dt/6*(dudt1+2*dudt2+2*dudt3+dudt4);
    v = v+dt/6*(dvdt1+2*dvdt2+2*dvdt3+dvdt4);
end
if nargin == 6
    err(end) = norm(ref(t(end))-u,'fro');
end

end
