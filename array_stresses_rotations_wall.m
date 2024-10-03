function [sigma,Omega]= array_stresses_rotations_wall(x1,x2,X2,P,A,CS,Q0,h)

epsilon=zeros(3,3);
u_g=zeros(3,3);
w=zeros(3,3);

M_P=diag(P);

%Gradient of displacement   
g_sigma=g0_wall(x1,x2,X2,Q0,P,P,h);
    
u_g(:,1)=A*g_sigma;
u_g(:,2)=A*(M_P*g_sigma);
u_g(:,3)=[0;0;0];
%Strain+rotation
for i=1:3
    for j=1:3
        epsilon(i,j)=1/2*(u_g(i,j)+u_g(j,i));
        w(i,j)=1/2*(u_g(i,j)-u_g(j,i));
    end
end
sigma=2*real(CS*[epsilon(1,1);epsilon(2,2);epsilon(3,3);2*epsilon(2,3);2*epsilon(3,1);2*epsilon(1,2)]);
Omega=2*real([-w(2,3);w(1,3);-w(1,2)]);
     







    
    