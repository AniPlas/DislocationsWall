function [sigma,Omega,strain]= array_stresses_rotations_bi_wall(x1,x2,X2,...
                P_I,P_conj_I,P_II,P_conj_II,A_I,A_II,CS_I,CS_II,Q01,Q02,W_I_II,V_conj_I_II,W_II_I,V_conj_II_I,h)

epsilon=zeros(3,3);
u_g=zeros(3,3);
w=zeros(3,3);

g1_sigma=zeros(3,1);
g2_sigma=zeros(3,1);

M_P_I=diag(P_I);
M_P_II=diag(P_II);

if (X2>=0) %Dislocations in Media I
    if (x2>0) %Point in Media I
    %Gradient of displacement
        g1_sigma(1)=V_conj_I_II(1,:)*conj(g0_wall(x1,x2,X2,Q01,P_conj_I(1),P_I,h));
        g1_sigma(2)=V_conj_I_II(2,:)*conj(g0_wall(x1,x2,X2,Q01,P_conj_I(2),P_I,h));
        g1_sigma(3)=V_conj_I_II(3,:)*conj(g0_wall(x1,x2,X2,Q01,P_conj_I(3),P_I,h));

        g_sigma=g0_wall(x1,x2,X2,Q01,P_I,P_I,h)+g1_sigma;

        u_g(:,1)=A_I*g_sigma;                                  
        u_g(:,2)=A_I*(M_P_I*g_sigma);  
        u_g(:,3)=[0;0;0]; 
    %Strain+rotation
        for i=1:3
            for j=1:3
                epsilon(i,j)=1/2*(u_g(i,j)+u_g(j,i));
                w(i,j)=1/2*(u_g(i,j)-u_g(j,i));
            end
        end
        sigma=2*real(CS_I*[epsilon(1,1);epsilon(2,2);epsilon(3,3);2*epsilon(2,3);2*epsilon(3,1);2*epsilon(1,2)]);
        Omega=2*real([-w(2,3);w(1,3);-w(1,2)]);
        strain=2*real(epsilon);

    else %Point in Media II
    %Gradient of displacement
        g2_sigma(1)=W_I_II(1,:)*g0_wall(x1,x2,X2,Q01,P_II(1),P_I,h);
        g2_sigma(2)=W_I_II(2,:)*g0_wall(x1,x2,X2,Q01,P_II(2),P_I,h); 
        g2_sigma(3)=W_I_II(3,:)*g0_wall(x1,x2,X2,Q01,P_II(3),P_I,h); 
        g_sigma=g2_sigma;    
        u_g(:,1)=A_II*g_sigma;                                  
        u_g(:,2)=A_II*(M_P_II*g_sigma);  
    %Strain+rotation
        for i=1:3
            for j=1:3
                epsilon(i,j)=1/2*(u_g(i,j)+u_g(j,i));
                w(i,j)=1/2*(u_g(i,j)-u_g(j,i));
            end
        end
        sigma=2*real(CS_II*[epsilon(1,1);epsilon(2,2);epsilon(3,3);2*epsilon(2,3);2*epsilon(3,1);2*epsilon(1,2)]);
        Omega=2*real([-w(2,3);w(1,3);-w(1,2)]);
        strain=2*real(epsilon);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %Dislocations in Media II
    if (x2<0) %Point in Media II
    %Gradient of displacement
        g1_sigma(1)=V_conj_II_I(1,:)*conj(g0_wall(x1,x2,X2,Q02,P_conj_II(1),P_II,h));
        g1_sigma(2)=V_conj_II_I(2,:)*conj(g0_wall(x1,x2,X2,Q02,P_conj_II(2),P_II,h));
        g1_sigma(3)=V_conj_II_I(3,:)*conj(g0_wall(x1,x2,X2,Q02,P_conj_II(3),P_II,h));

        g_sigma=g0_wall(x1,x2,X2,Q02,P_II,P_II,h)+g1_sigma;

        u_g(:,1)=A_II*g_sigma;                                  
        u_g(:,2)=A_II*(M_P_II*g_sigma);  
        u_g(:,3)=[0;0;0]; 
    %Strain+rotation
        for i=1:3
            for j=1:3
                epsilon(i,j)=1/2*(u_g(i,j)+u_g(j,i));
                w(i,j)=1/2*(u_g(i,j)-u_g(j,i));
            end
        end
        sigma=2*real(CS_II*[epsilon(1,1);epsilon(2,2);epsilon(3,3);2*epsilon(2,3);2*epsilon(3,1);2*epsilon(1,2)]);
        Omega=2*real([-w(2,3);w(1,3);-w(1,2)]);
        strain=2*real(epsilon);

    else %Point in Media I
    %Gradient of displacement
        g2_sigma(1)=W_II_I(1,:)*g0_wall(x1,x2,X2,Q02,P_I(1),P_II,h);
        g2_sigma(2)=W_II_I(2,:)*g0_wall(x1,x2,X2,Q02,P_I(2),P_II,h); 
        g2_sigma(3)=W_II_I(3,:)*g0_wall(x1,x2,X2,Q02,P_I(3),P_II,h); 
        g_sigma=g2_sigma;    
        u_g(:,1)=A_I*g_sigma;                                  
        u_g(:,2)=A_I*(M_P_I*g_sigma);  
    %Strain+rotation
        for i=1:3
            for j=1:3
                epsilon(i,j)=1/2*(u_g(i,j)+u_g(j,i));
                w(i,j)=1/2*(u_g(i,j)-u_g(j,i));
            end
        end
        sigma=2*real(CS_I*[epsilon(1,1);epsilon(2,2);epsilon(3,3);2*epsilon(2,3);2*epsilon(3,1);2*epsilon(1,2)]);
        Omega=2*real([-w(2,3);w(1,3);-w(1,2)]);
        strain=2*real(epsilon);
    end
end







    
    