% Author: Thiebaud RICHETON
% Date: October 2024
% email: thiebaud.richeton@univ-lorraine.fr

clear;
close all;
clc;

% Material='Cu';
c11=169.886;  % in GPa
c12=122.609;  % in GPa
c44=76.191;  % in GPa
a=0.3615; % cell parameters in nm

Azener=2*c44/(c11-c12) % Anisotropy ratio 

% Elastic constants for comparisons with isotropic elasticity
% Computation of Voigt-Reuss-Hill average
mu1=0.5*(c11-c12);
mu2=c44;
GV=2/5*mu1+3/5*mu2; % Voigt model for shear modulus
GR=5/(2/mu1+3/mu2); % Reuss model for shear modulus
mu = 0.5*(GV+GR) % Voigt-Reuss-Hill average
K=1/3*c11+2/3*c12; % Voigt and Reuss models for bulk modulus
nu = (3*K-2*mu)/(2*(3*K+mu)) % Poisson's ratio
E = 2*mu*(1+nu) % Young modulus

b=sqrt(2)/2*a % perfect Burgers vector magnitude for FCC in nm

h=2.81 % Distance between dislocation in the peridoci GB wall (in nm)

x = [1;0;0];
y = [0;1;0];
z = [0;0;1];
I=diag([1,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bo=b/sqrt(2)*[-1;1;0]; % Burgers vector in crystal coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge dislocations
to=[1;1;2]; % Dislocation line direction in crystal coordinates
% to=[0;0;1]; % Dislocation line direction in crystal coordinates
to=to/norm(to);
nGBo=bo; % GB boundary normal in crystal coordinates
nGBo=nGBo/norm(nGBo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Screw dislocations
% to=bo; % Dislocation line direction in crystal coordinates
% to=to/norm(to);
% nGBo=[1;1;2]; % GB boundary normal in crystal coordinates
% nGBo=nGBo/norm(nGBo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zo=cross(nGBo,to);
zo=zo/norm(zo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=[0;0;1]; % Dislocation line direction vector in global coordinates
nGB0=[0;1;0]; % GB boundary normal in global coordinates
z0=cross(nGB0,t0);
z0=z0/norm(z0);

% Transition matrix from crystal coordinates to global coordinates
T=[nGB0,t0,z0]/[nGBo,to,zo];
    
bGB=T*bo; % Burgers vector in global coordinates

% Construct elastic constant, C is in form tetradic, CS is in form 6*6 matrix (the contracted Voigt notation)
[C,CS]=elasticconstcubic(c11,c12,c44,T);
% Construct matrix P, A and B using Stroh matrix Formalism, _conj means conjugate of the matrix, M_P_=diag(P_), M_=i*A_/B_
[P,P_conj,A,A_conj,B,M,M_conj]=Choi(C);
Q0=-0.5*(-1)^0.5/pi*(B.'*bGB);
M_P=diag(P);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of computation in nm for 1st plot check
x_1=0.5*h;
dd=10*h;
x_2=linspace(-dd,dd,1E3);
profile=x_2;

Rotation1=zeros(length(profile),1);
Rotation2=zeros(length(profile),1);
Rotation3=zeros(length(profile),1);
Omega1walliso=zeros(length(profile),1);
Omega2walliso=zeros(length(profile),1);
Omega3walliso=zeros(length(profile),1);
% Compute stresses and rotations
for i=1:length(profile)
    Z=P*x_2(i)+x_1;
    [sigma,Omega]=array_stresses_rotations_wall(x_1,x_2(i),0,P,A,CS,Q0,h);
    Rotation1(i)=Omega(1);
    Rotation2(i)=Omega(2);
    Rotation3(i)=Omega(3);
    % Edge dislocations
    Omega3walliso(i)=-b/(2*h)*sinh(2*pi*x_2(i)/h)/(cosh(2*pi*x_2(i)/h)-cos(2*pi*x_1/h));
    % Screw dislocations
%     Omega1walliso(i)=b/(4*h)*sin(2*pi*x_1/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1/h));
%     Omega2walliso(i)=b/(4*h)*sinh(2*pi*(x_2(i))/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1/h));
    if (x_2(i)<-9*h)
        kII=i;
    end
    if (x_2(i)<9*h)
        kI=i+1;
    end
end

Omega1_I=mean(Rotation1(kI:length(x_2)));
Omega1_II=mean(Rotation1(1:kII));
Omega2_I=mean(Rotation2(kI:length(x_2)));
Omega2_II=mean(Rotation2(1:kII));
Omega3_I=mean(Rotation3(kI:length(x_2)));
Omega3_II=mean(Rotation3(1:kII));

omega_I=[[0 -Omega3_I Omega2_I];...
         [Omega3_I 0 -Omega1_I];...
         [-Omega2_I Omega1_I 0]];
     
omega_II=[[0 -Omega3_II Omega2_II];...
         [Omega3_II 0 -Omega1_II];...
         [-Omega2_II Omega1_II 0]];

% Edge dislocations
Omega3iso_I=mean(Omega3walliso(kI:length(x_2)));
Omega3iso_II=mean(Omega3walliso(1:kII)); 
% screw dislocations
% Omega2iso_I=mean(Omega2walliso(kI:length(x_2)));
% Omega2iso_II=mean(Omega2walliso(1:kII)); 
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition matrix from first crystal coordinates to global coordinates
rx=omega_I(2,3)-omega_I(3,2);
ry=omega_I(3,1)-omega_I(1,3);
rz=omega_I(1,2)-omega_I(2,1);
N=sqrt(rx^2+ry^2+rz^2);
axe=[rx/N;ry/N;rz/N]
theta=asin(0.5*N);
theta1_deg=theta*180/pi
R_I=mrot(axe,theta);

% Transition matrix from second crystal coordinates to global coordinates
rx=omega_II(2,3)-omega_II(3,2);
ry=omega_II(3,1)-omega_II(1,3);
rz=omega_II(1,2)-omega_II(2,1);
N=sqrt(rx^2+ry^2+rz^2);
axe=[rx/N;ry/N;rz/N]
theta=asin(0.5*N);
theta2_deg=theta*180/pi
R_II=mrot(axe,theta);

theta_aniso=theta1_deg+theta2_deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
set(gca,'FontSize',20)
plot(x_2,Rotation1.*180/pi,'-b','LineWidth',3)
plot(x_2,Omega1walliso.*180/pi,':b','LineWidth',3)
plot(x_2,Rotation2.*180/pi,'-g','LineWidth',3)
plot(x_2,Omega2walliso.*180/pi,':g','LineWidth',3)
plot(x_2,Rotation3.*180/pi,'-r','LineWidth',3)
plot(x_2,Omega3walliso.*180/pi,':r','LineWidth',3)
xlabel('profile (nm)')
ylabel('\Omega (°)')
legend('\Omega_1','\Omega_1 iso','\Omega_2','\Omega_2 iso','\Omega_3','\Omega_3 iso')
legend('boxoff')
title('check at first step')
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of computation in nm for rotation convergence
x_1=0.5*h;
dd=10*h;
A=linspace(-10*h,-9*h,1E1);
B=linspace(9*h,10*h,1E1);
x_2=[A B];
profile=x_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_aniso_old=theta_aniso;
theta_aniso_new=2*theta_aniso;
while (abs(theta_aniso_new-theta_aniso_old)/theta_aniso_new*100>1E-3)
    theta_aniso_old=theta_aniso;
% Construct elastic constant, C_ is in form tetradic, CS_ is in form 6*6 matrix (the contracted Voigt notation)
    [C_I,CS_I]=elasticconstmis(C,R_I);
    [C_II,CS_II]=elasticconstmis(C,R_II);

% Construct matrix P, A and B using Stroh matrix Formalism, _conj means conjugate of the matrix, M_P_=diag(P_), M_=i*A_/B_
    [P_I,P_conj_I,A_I,A_conj_I,B_I,M_I,M_conj_I]=Choi(C_I);
    [P_II,P_conj_II,A_II,A_conj_II,B_II,M_II,M_conj_II]=Choi(C_II);

    T_I_II=(M_II+conj(M_I))\(M_I-M_II);
    W_I_II=B_II\(I+T_I_II)*B_I;
    V_I_II=conj(B_I)\T_I_II*B_I;
    V_conj_I_II=conj(V_I_II);

    T_II_I=(M_I+conj(M_II))\(M_II-M_I);
    W_II_I=B_I\(I+T_II_I)*B_II;
    V_II_I=conj(B_II)\T_II_I*B_II;
    V_conj_II_I=conj(V_II_I);

    Q0=-0.5*(-1)^0.5/pi*(B_I.'*bGB);

% Compute stresses and rotations
    for i=1:length(profile)
        [sigma,Omega,strain]=array_stresses_rotations_bi_wall(x_1,x_2(i),0,...
            P_I,P_conj_I,P_II,P_conj_II,A_I,A_II,CS_I,CS_II,Q0,Q0,W_I_II,V_conj_I_II,W_II_I,V_conj_II_I,h);
        Rotation1(i)=Omega(1);
        Rotation2(i)=Omega(2);
        Rotation3(i)=Omega(3);

        if (x_2(i)<-9*h)
            kII=i;
        end
        if (x_2(i)<9*h)
            kI=i+1;
        end
    end

    Omega1_I=mean(Rotation1(kI:length(x_2)));
    Omega1_II=mean(Rotation1(1:kII));
    Omega2_I=mean(Rotation2(kI:length(x_2)));
    Omega2_II=mean(Rotation2(1:kII));
    Omega3_I=mean(Rotation3(kI:length(x_2)));
    Omega3_II=mean(Rotation3(1:kII));

    omega_I=[[0 -Omega3_I Omega2_I];...
             [Omega3_I 0 -Omega1_I];...
             [-Omega2_I Omega1_I 0]];

    omega_II=[[0 -Omega3_II Omega2_II];...
             [Omega3_II 0 -Omega1_II];...
             [-Omega2_II Omega1_II 0]];
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition matrix from first crystal coordinates to global coordinates
    rx=omega_I(2,3)-omega_I(3,2);
    ry=omega_I(3,1)-omega_I(1,3);
    rz=omega_I(1,2)-omega_I(2,1);
    N=sqrt(rx^2+ry^2+rz^2);
    axe=[rx/N;ry/N;rz/N];
    theta=asin(0.5*N);
    theta1_deg=theta*180/pi;
    R_I=mrot(axe,theta);

% Transition matrix from second crystal coordinates to global coordinates
    rx=omega_II(2,3)-omega_II(3,2);
    ry=omega_II(3,1)-omega_II(1,3);
    rz=omega_II(1,2)-omega_II(2,1);
    N=sqrt(rx^2+ry^2+rz^2);
    axe=[rx/N;ry/N;rz/N];
    theta=asin(0.5*N);
    theta2_deg=theta*180/pi;
    R_II=mrot(axe,theta);

    theta_aniso=theta1_deg+theta2_deg
    theta_aniso_new=theta_aniso;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of computation in nm for profile plot
% x_1=linspace(0.1*h,0.9*h,1E2);
% x_2=0*h*ones(length(x_1),1);
% profile=x_1;

x_2=linspace(-10*h,10*h,1E3);
x_1=0.8*h*ones(length(x_2),1);
profile=x_2;

Rotation1=zeros(length(profile),1);
Rotation2=zeros(length(profile),1);
Rotation3=zeros(length(profile),1);
Sig11=zeros(length(profile),1);
Sig22=zeros(length(profile),1);
Sig33=zeros(length(profile),1);
Sig23=zeros(length(profile),1);
Sig31=zeros(length(profile),1);
Sig12=zeros(length(profile),1);
Strain11=zeros(length(profile),1);
Strain22=zeros(length(profile),1);
Strain33=zeros(length(profile),1);
Strain23=zeros(length(profile),1);
Strain31=zeros(length(profile),1);
Strain12=zeros(length(profile),1);
we=zeros(length(profile),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct elastic constant, C_ is in form tetradic, CS_ is in form 6*6 matrix (the contracted Voigt notation)
[C_I,CS_I]=elasticconstmis(C,R_I);
[C_II,CS_II]=elasticconstmis(C,R_II);

% Construct matrix P, A and B using Stroh matrix Formalism, _conj means conjugate of the matrix, M_P_=diag(P_), M_=i*A_/B_
[P_I,P_conj_I,A_I,A_conj_I,B_I,M_I,M_conj_I]=Choi(C_I);
[P_II,P_conj_II,A_II,A_conj_II,B_II,M_II,M_conj_II]=Choi(C_II);

T_I_II=(M_II+conj(M_I))\(M_I-M_II);
W_I_II=B_II\(I+T_I_II)*B_I;
V_I_II=conj(B_I)\T_I_II*B_I;
V_conj_I_II=conj(V_I_II);
M_P_I=diag(P_I);

T_II_I=(M_I+conj(M_II))\(M_II-M_I);
W_II_I=B_I\(I+T_II_I)*B_II;
V_II_I=conj(B_II)\T_II_I*B_II;
V_conj_II_I=conj(V_II_I);

Q0=-0.5*(-1)^0.5/pi*(B_I.'*bGB);

% Compute stresses and rotations
for i=1:length(profile)
    [sigma,Omega,strain]=array_stresses_rotations_bi_wall(x_1(i),x_2(i),0,...
        P_I,P_conj_I,P_II,P_conj_II,A_I,A_II,CS_I,CS_II,Q0,Q0,W_I_II,V_conj_I_II,W_II_I,V_conj_II_I,h);
    Rotation1(i)=Omega(1);
    Rotation2(i)=Omega(2);
    Rotation3(i)=Omega(3);
    Sig11(i)=sigma(1);
    Sig22(i)=sigma(2);
    Sig33(i)=sigma(3);
    Sig23(i)=sigma(4);
    Sig31(i)=sigma(5);
    Sig12(i)=sigma(6);
    Strain11(i)=strain(1,1);
    Strain22(i)=strain(2,2);
    Strain33(i)=strain(3,3);
    Strain23(i)=strain(2,3);
    Strain31(i)=strain(3,1);
    Strain12(i)=strain(1,2);
    
    we(i)=0.5*(Sig11(i)*Strain11(i)+Sig22(i)*Strain22(i)+Sig33(i)*Strain33(i)+...
            2*Sig23(i)*Strain23(i)+2*Sig31(i)*Strain31(i)+2*Sig12(i)*Strain12(i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute infinite wall stresses and rotations in isotropic elasticity and stresses for finite wall in anisotropic elasticity
Sig_iso_wall_11=zeros(length(profile),1);
Sig_iso_wall_22=zeros(length(profile),1);
Sig_iso_wall_33=zeros(length(profile),1);
Sig_iso_wall_12=zeros(length(profile),1);
Sig_iso_wall_23=zeros(length(profile),1);
Sig_iso_wall_31=zeros(length(profile),1);
Omega1walliso=zeros(length(profile),1);
Omega2walliso=zeros(length(profile),1);
Omega3walliso=zeros(length(profile),1);
wewalliso=zeros(length(profile),1);
for i=1:length(profile)
    ksi=pi*x_2(i)/h;
    lambda=pi*x_1(i)/h;
    
    % Edge dislocations
    D=(mu*b)/(2*h*(1-nu)*(cosh(2*ksi)-cos(2*lambda))^2);
    
    Sig_iso_wall_11(i)=((sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)-2*ksi*sinh(2*ksi)))*D;
    Sig_iso_wall_22(i)=((sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)+2*ksi*sinh(2*ksi)))*D;
    Sig_iso_wall_33(i)=((2*nu*sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)))*D;
    Sig_iso_wall_12(i)=-((2*ksi)*(cosh(2*ksi)*cos(2*lambda)-1))*D;
    
    % Screw dislocations
%     Sig_iso_wall_23(i)=mu*b/(2*h)*sin(2*lambda)/(cosh(2*ksi)-cos(2*lambda));
%     Sig_iso_wall_31(i)=-mu*b/(2*h)*sinh(2*ksi)/(cosh(2*ksi)-cos(2*lambda)); 
    
    % Edge dislocations
    Omega3walliso(i)=-b/(2*h)*sinh(2*pi*x_2(i)/h)/(cosh(2*pi*x_2(i)/h)-cos(2*pi*x_1(i)/h));
    % Screw dislocations
%     Omega1walliso(i)=b/(4*h)*sin(2*pi*x_1(i)/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1(i)/h));
%     Omega2walliso(i)=b/(4*h)*sinh(2*pi*(x_2(i))/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1(i)/h));

    Strain_iso_wall_11=1/E*(Sig_iso_wall_11(i)-nu*(Sig_iso_wall_22(i)+Sig_iso_wall_33(i)));
    Strain_iso_wall_22=1/E*(Sig_iso_wall_22(i)-nu*(Sig_iso_wall_11(i)+Sig_iso_wall_33(i)));
    Strain_iso_wall_33=1/E*(Sig_iso_wall_33(i)-nu*(Sig_iso_wall_11(i)+Sig_iso_wall_22(i)));
    Strain_iso_wall_23=1/(2*mu)*Sig_iso_wall_23(i);
    Strain_iso_wall_31=1/(2*mu)*Sig_iso_wall_31(i);
    Strain_iso_wall_12=1/(2*mu)*Sig_iso_wall_12(i);
    
    wewalliso(i)=0.5*(Sig_iso_wall_11(i)*Strain_iso_wall_11+Sig_iso_wall_22(i)*Strain_iso_wall_22+...
        Sig_iso_wall_33(i)*Strain_iso_wall_33+2*Sig_iso_wall_23(i)*Strain_iso_wall_23+...
        2*Sig_iso_wall_31(i)*Strain_iso_wall_31+2*Sig_iso_wall_12(i)*Strain_iso_wall_12);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
hold on
set(gca,'FontSize',20)
plot(profile,Sig11,'-b','LineWidth',3)
plot(profile,Sig_iso_wall_11,':b','LineWidth',3)
plot(profile,Sig22,'-g','LineWidth',3)
plot(profile,Sig_iso_wall_22,':g','LineWidth',3)
plot(profile,Sig33,'-r','LineWidth',3)
plot(profile,Sig_iso_wall_33,':r','LineWidth',3)
plot(profile,Sig11+Sig22+Sig33,'-k','LineWidth',3)
plot(profile,Sig_iso_wall_11+Sig_iso_wall_22+Sig_iso_wall_33,':k','LineWidth',3)
xlabel('profile (nm)')
ylabel('\sigma (GPa)')
legend('\sigma_{11}','\sigma_{11} iso','\sigma_{22}','\sigma_{22} iso','\sigma_{33}','\sigma_{33} iso','\sigma_{kk}','\sigma_{kk} iso')
legend('boxoff')
box on

figure(3)
hold on
set(gca,'FontSize',20)
plot(profile,Sig12,'-b','LineWidth',3)
plot(profile,Sig_iso_wall_12,':b','LineWidth',3)
plot(profile,Sig31,'-g','LineWidth',3)
plot(profile,Sig_iso_wall_31,':g','LineWidth',3)
plot(profile,Sig23,'-r','LineWidth',3)
plot(profile,Sig_iso_wall_23,':r','LineWidth',3)
xlabel('profile (nm)')
ylabel('\sigma (GPa)')
legend('\sigma_{12}','\sigma_{12} iso','\sigma_{31}','\sigma_{31} iso','\sigma_{23}','\sigma_{23} iso')
legend('boxoff')
box on

figure(10)
hold on
set(gca,'FontSize',20)
plot(profile,Rotation1.*180/pi,'-b','LineWidth',3)
plot(profile,Omega1walliso.*180/pi,':b','LineWidth',3)
plot(profile,Rotation2.*180/pi,'-g','LineWidth',3)
plot(profile,Omega2walliso.*180/pi,':g','LineWidth',3)
plot(profile,Rotation3.*180/pi,'-r','LineWidth',3)
plot(profile,Omega3walliso.*180/pi,':r','LineWidth',3)
xlabel('profile (nm)')
ylabel('\Omega (°)')
legend('\Omega_1','\Omega_1 iso','\Omega_2','\Omega_2 iso','\Omega_3','\Omega_3 iso')
legend('boxoff')
box on

figure(7)
hold on
set(gca,'FontSize',20)
plot(profile,Strain11,'b','LineWidth',3)
plot(profile,Strain22,'g','LineWidth',3)
plot(profile,Strain33,'r','LineWidth',3)
plot(profile,Strain23,'k','LineWidth',3)
plot(profile,Strain31,'m','LineWidth',3)
plot(profile,Strain12,'y','LineWidth',3)
xlabel('profile (nm)')
ylabel('Strain')
legend('\epsilon_{11}','\epsilon_{22}','\epsilon_{33}','\epsilon_{23}','\epsilon_{31}','\epsilon_{12}')
legend('boxoff')
box on

figure(8)
hold on
set(gca,'FontSize',20)
plot(profile,we,'-r','LineWidth',3)
plot(profile,wewalliso,'--b','LineWidth',3)
xlabel('profile (nm)')
ylabel('\omega_e')
legend('aniso','iso')
legend('boxoff')
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positions of computation in nm for contour plots
x_1=linspace(0*h,1*h,1E2);
x_2=linspace(-1*h,1*h,1E2);

Rotation1=zeros(length(x_2),length(x_1));
Rotation2=zeros(length(x_2),length(x_1));
Rotation3=zeros(length(x_2),length(x_1));
Sig11=zeros(length(x_2),length(x_1));
Sig22=zeros(length(x_2),length(x_1));
Sig33=zeros(length(x_2),length(x_1));
Sig23=zeros(length(x_2),length(x_1));
Sig31=zeros(length(x_2),length(x_1));
Sig12=zeros(length(x_2),length(x_1));
Strain11=zeros(length(x_2),length(x_1));
Strain22=zeros(length(x_2),length(x_1));
Strain33=zeros(length(x_2),length(x_1));
Strain23=zeros(length(x_2),length(x_1));
Strain31=zeros(length(x_2),length(x_1));
Strain12=zeros(length(x_2),length(x_1));
we=zeros(length(x_2),length(x_1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stresses and rotations
for i=1:length(x_2)
    for j=1:length(x_1)
        [sigma,Omega,strain]=array_stresses_rotations_bi_wall(x_1(j),x_2(i),0,...
            P_I,P_conj_I,P_II,P_conj_II,A_I,A_II,CS_I,CS_II,Q0,Q0,W_I_II,V_conj_I_II,W_II_I,V_conj_II_I,h);
        Rotation1(j,i)=Omega(1);
        Rotation2(j,i)=Omega(2);
        Rotation3(j,i)=Omega(3);
        Sig11(j,i)=sigma(1);
        Sig22(j,i)=sigma(2);
        Sig33(j,i)=sigma(3);
        Sig23(j,i)=sigma(4);
        Sig31(j,i)=sigma(5);
        Sig12(j,i)=sigma(6);
        Strain11(j,i)=strain(1,1);
        Strain22(j,i)=strain(2,2);
        Strain33(j,i)=strain(3,3);
        Strain23(j,i)=strain(2,3);
        Strain31(j,i)=strain(3,1);
        Strain12(j,i)=strain(1,2);
    
        we(j,i)=0.5*(Sig11(j,i)*Strain11(j,i)+Sig22(j,i)*Strain22(j,i)+Sig33(j,i)*Strain33(j,i)+...
            2*Sig23(j,i)*Strain23(j,i)+2*Sig31(j,i)*Strain31(j,i)+2*Sig12(j,i)*Strain12(j,i));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute infinite wall stresses and rotations in isotropic elasticity and stresses for finite wall in anisotropic elasticity
Sig_iso_wall_11=zeros(length(x_2),length(x_1));
Sig_iso_wall_22=zeros(length(x_2),length(x_1));
Sig_iso_wall_33=zeros(length(x_2),length(x_1));
Sig_iso_wall_12=zeros(length(x_2),length(x_1));
Sig_iso_wall_23=zeros(length(x_2),length(x_1));
Sig_iso_wall_31=zeros(length(x_2),length(x_1));
Omega1walliso=zeros(length(x_2),length(x_1));
Omega2walliso=zeros(length(x_2),length(x_1));
Omega3walliso=zeros(length(x_2),length(x_1));
wewalliso=zeros(length(x_2),length(x_1));

for i=1:length(x_2)
    for j=1:length(x_1)
        ksi=pi*x_2(i)/h;
        lambda=pi*x_1(j)/h;
    
% Edge dislocations
        D=(mu*b)/(2*h*(1-nu)*(cosh(2*ksi)-cos(2*lambda))^2);

        Sig_iso_wall_11(j,i)=((sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)-2*ksi*sinh(2*ksi)))*D;
        Sig_iso_wall_22(j,i)=((sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)+2*ksi*sinh(2*ksi)))*D;
        Sig_iso_wall_33(j,i)=((2*nu*sin(2*lambda))*(cosh(2*ksi)-cos(2*lambda)))*D;
        Sig_iso_wall_12(j,i)=-((2*ksi)*(cosh(2*ksi)*cos(2*lambda)-1))*D;
    
% Screw dislocations
    Sig_iso_wall_23(j,i)=mu*b/(2*h)*sin(2*lambda)/(cosh(2*ksi)-cos(2*lambda));
    Sig_iso_wall_31(j,i)=-mu*b/(2*h)*sinh(2*ksi)/(cosh(2*ksi)-cos(2*lambda)); 
    
% Edge dislocations
        Omega3walliso(j,i)=-b/(2*h)*sinh(2*pi*x_2(i)/h)/(cosh(2*pi*x_2(i)/h)-cos(2*pi*x_1(j)/h));
% Screw dislocations
%     Omega1walliso(j,i)=b/(4*h)*sin(2*pi*x_1(j)/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1(j)/h));
%     Omega2walliso(j,i)=b/(4*h)*sinh(2*pi*(x_2(i))/h)/(cosh(2*pi*(x_2(i))/h)-cos(2*pi*x_1(j)/h));

    Strain_iso_wall_11=1/E*(Sig_iso_wall_11(j,i)-nu*(Sig_iso_wall_22(j,i)+Sig_iso_wall_33(j,i)));
    Strain_iso_wall_22=1/E*(Sig_iso_wall_22(j,i)-nu*(Sig_iso_wall_11(j,i)+Sig_iso_wall_33(j,i)));
    Strain_iso_wall_33=1/E*(Sig_iso_wall_33(j,i)-nu*(Sig_iso_wall_11(j,i)+Sig_iso_wall_22(j,i)));
    Strain_iso_wall_23=1/(2*mu)*Sig_iso_wall_23(j,i);
    Strain_iso_wall_31=1/(2*mu)*Sig_iso_wall_31(j,i);
    Strain_iso_wall_12=1/(2*mu)*Sig_iso_wall_12(j,i);

    wewalliso(j,i)=0.5*(Sig_iso_wall_11(j,i)*Strain_iso_wall_11+Sig_iso_wall_22(j,i)*Strain_iso_wall_22+...
        Sig_iso_wall_33(j,i)*Strain_iso_wall_33+2*Sig_iso_wall_23(j,i)*Strain_iso_wall_23+...
        2*Sig_iso_wall_31(j,i)*Strain_iso_wall_31+2*Sig_iso_wall_12(j,i)*Strain_iso_wall_12);
    end
end

V_legend=0.5*(-10:1:10);
figure(100)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig11,V_legend,'Showtext','on');
colorbar;
title("\sigma_{11} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(101)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig22,V_legend,'Showtext','on');
colorbar;
title("\sigma_{22} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(102)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig33,V_legend,'Showtext','on');
colorbar;
title("\sigma_{33} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(103)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig23,V_legend,'Showtext','on');
colorbar;
title("\sigma_{23} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(104)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig31,V_legend,'Showtext','on');
colorbar;
title("\sigma_{31} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(105)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig12,V_legend,'Showtext','on');
colorbar;
title("\sigma_{12} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(110)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig11+Sig22+Sig33,V_legend,'Showtext','on');
colorbar;
title("\sigma_{kk} (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(120)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,we,0.2*(0:0.1:1),'Showtext','on');
colorbar;
title("\omega_e (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(200)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_11,V_legend,'Showtext','on');
colorbar;
title("\sigma_{11} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(201)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_22,V_legend,'Showtext','on');
colorbar;
title("\sigma_{22} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(202)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_33,V_legend,'Showtext','on');
colorbar;
title("\sigma_{33} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(203)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_23,V_legend,'Showtext','on');
colorbar;
title("\sigma_{23} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(204)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_31,V_legend,'Showtext','on');
colorbar;
title("\sigma_{31} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(205)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_12,V_legend,'Showtext','on');
colorbar;
title("\sigma_{12} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(210)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,Sig_iso_wall_11+Sig_iso_wall_22+Sig_iso_wall_33,V_legend,'Showtext','on');
colorbar;
title("\sigma_{kk} iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on

figure(220)
set(gca,'FontSize',20)
contourf(x_2/b,x_1/b,wewalliso,0.2*(0:0.1:1),'Showtext','on');
colorbar;
title("\omega_e iso (GPa)");
xlabel("x_2");
ylabel("x_1");
box on




