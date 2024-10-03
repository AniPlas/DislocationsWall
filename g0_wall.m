function g0_wall=g0_wall(x1,x2,X2,Q0,Pz,PS,h)
g0_wall=-pi/h*diag(cot((PS*X2-Pz*x2-x1)*pi/h))*Q0;
