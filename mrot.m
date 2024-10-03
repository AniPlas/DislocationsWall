function [R] = mrot(X,theta)

R = [ ];

s=theta;

x=X(1);
y=X(2);
z=X(3);

%%%%%%%%%%%%%%%%%%%%%%%%%% rotation clockwise %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% angles Euler EBSD : passage du repere global au repere du cristal par des
% rotations anti-horaires (convention Bunge)
% --> pour determiner la matrice de passage global --> cristal, on ecrit
% les matrices de rotation dans le sens horaire

R(1,1)=cos(s)+x^2*(1-cos(s));

R(1,2)=x*y*(1-cos(s))+z*sin(s);

R(1,3)=x*z*(1-cos(s))-y*sin(s);

R(2,1)=y*x*(1-cos(s))-z*sin(s);

R(2,2)=cos(s)+y^2*(1-cos(s));

R(2,3)=y*z*(1-cos(s))+x*sin(s);

R(3,1)=z*x*(1-cos(s))+y*sin(s);

R(3,2)=z*y*(1-cos(s))-x*sin(s);

R(3,3)=cos(s)+z^2*(1-cos(s));


