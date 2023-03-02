function [z , v] = BS (alpha,beta,gamma,n,r)
%% Backward Substitution For Solvig Subsystems A11 Z = p and A11 z = b_2
alpha_1=alpha/beta;
gamma_1=gamma/beta;
z=zeros(1,n-1);
z(n-1)=r(n-1)/beta;
z(n-2)=r(n-2)/beta-alpha_1*z(n-1);
for i=2:n-2
    z(n-i-1)=(r(n-3)/beta)-alpha_1*z(n-i)-gamma_1*z(n-i+1);
end


