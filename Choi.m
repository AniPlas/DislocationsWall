function [P,P_conj,A,A_conj,B,M,M_conj]=Choi(C)
% calculate P, A and B using Stroh matrix Formalism
% calculate QI and M by equation (6) and (8)
% C: (3x3x3x3, rank-4 tensor) elastic constant tensor

%Define Q, R and T matrices 
z3 = zeros(3);Q=z3; R=z3; T=z3;
for i=1:3
     for k=1:3
                 Q(i,k)=C(i,1,k,1);
                 R(i,k)=C(i,1,k,2);
                 T(i,k)=C(i,2,k,2);
     end
end
 
%Create 6x6 matrix N0
%By theory of T. C. T. Ting. 
%ACTA MECHANICA SINICA, Vol.8, No.3, 193-207, August 1992. Generalized
%Stroh formalism for anisotropic elasticity for general boundary conditions
N0=zeros(6);
invT = inv(T);
N0=[-invT*R',invT;R/T*R'-Q,-R/T];
% N0=[-invT*R',invT;R*invT*R'-Q,-R*(invT)'];  %Ting %correction S. Berbenni
 
%Solve Eigen equation
[V,D]=eig(N0);
 
%initializing A, B and P
P=zeros(1,3); 
A = z3; B = z3; 

%Define A, B and P
for alpha=1:3
    if imag(D(2*alpha,2*alpha))>0
        P(alpha)=D(2*alpha,2*alpha);
        A(:,alpha)=V(1:3,2*alpha);
        B(:,alpha)=V(4:6,2*alpha);
    else
        P(alpha)=D(2*alpha-1,2*alpha-1);
        A(:,alpha)=V(1:3,2*alpha-1);
        B(:,alpha)=V(4:6,2*alpha-1);
    end
end

%normalize A and B
for alpha=1:3
    %AialphaLibeta+AibetaLialpha=delta(alpha,beta)
    ralpha = sqrt(2*sum(A(:,alpha).*B(:,alpha)));
    A(:,alpha) = A(:,alpha)/ralpha;
    B(:,alpha) = B(:,alpha)/ralpha;    
end


%Construct QI and M matrices using A, B and b by equation (6) and (8)
i=sqrt(-1);
M=i*A/B;

P_conj=conj(P);
A_conj=conj(A);
M_conj=conj(M);
P=P.';

