function a = int_d1(p)
% 2p-th order central difference for first derivative

nu = sym(-p:p);
a = (-1).^(nu+1)*factorial(p)^2./nu./factorial(p-nu)./factorial(p+nu);
a(p+1) = 0;

end
