function J=mop3(theta)
%  J=mop3(theta)
%  J=[J1 J2]
%  theta=[theta_1 theta_2]
%  -pi<=theta_i<pi i=1,2

a1=0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
a2=1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
b1=0.5*sin(theta(1,1))-2*cos(theta(1,1))+sin(theta(1,2))-1.5*cos(theta(1,2));
b2=1.5*sin(theta(1,1))-cos(theta(1,1))+2*sin(theta(1,2))-0.5*cos(theta(1,2));

J1=(1+(a1-b1)^2+(a2-b2)^2);
J2=((theta(1,1)+3)^2+(theta(1,2)+1)^2);
J=[J1 J2];
