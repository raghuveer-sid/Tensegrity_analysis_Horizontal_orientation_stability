clear all
syms theta z x w rf beta k F h alpha

%Rotation Matrices
   
Rzalpha = [cos(alpha) sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1];
Rybeta = [cos(beta) 0 sin(-beta);0 1 0;sin(beta) 0 cos(beta)];

R =Rzalpha*Rybeta;

%End effector coordinates

o1 = [-h*rf ;0 ;rf];
o2 = [-h*rf ;rf*sin(2*pi/3);rf*cos(2*pi/3)];
o3 = [-h*rf ;rf*sin(2*2*pi/3);rf*cos(2*2*pi/3)];
c1 = R*o1;
c2 = R*o2;
c3 = R*o3;

%Base coordinates

b1= [h*rf; 0; rf];
b2 = [h*rf;rf*sin(2*pi/3);rf*cos(2*pi/3)];
b3 = [h*rf;rf*sin(2*2*pi/3);rf*cos(2*2*pi/3)];

%Length of tensegrity

l1 = sqrt((b1(1)-c1(1))^2+(b1(2)-c1(2))^2+(b1(3)-c1(3))^2);
l2 = sqrt((b2(1)-c2(1))^2+(b2(2)-c2(2))^2+(b2(3)-c2(3))^2);
l3 = sqrt((b3(1)-c3(1))^2+(b3(2)-c3(2))^2+(b3(3)-c3(3))^2);

%Force vector

r1 = c2-c1;
r2 = c3-c1;
r_1 = norm(r1);
r_2 = norm(r2);
f1 = w*r1/r_1;
f2 = w*r2/r_2;

%Moments

mc = cross(c2,f1)+cross(c3,f2);
ms = k*(cross(c2,f1)+cross(c3,f2));

%Potential Energy

uc = F*(l1+l2+l3);
us = 0.5*k*(l1^2+l2^2+l3^2);
ut = uc+us;
kalpha= diff(ut,alpha,2);

alpha_1 = 0;
beta_1 = 0;
k_1 = 0.2;
F_1 = 1;
h_1 = 0.5;
rf_1 = [1:2:20]';

new_kalpha = subs(kalpha, {beta,k,F,h,rf,alpha},{beta_1,k_1,F_1,h_1,rf_1,alpha_1});

kalpha_val= convert_to_num(new_kalpha);


G = kalpha;
figure(2)
plot(rf_1,kalpha_val);
title('Kalpha vs rf')
xlabel('rf')
ylabel('Kalpha')

function num_array = convert_to_num(sym_Data)
num_array=[];
for i = 1:length(sym_Data)
    sym_Data_cell=sym2cell(sym_Data);
   num_array(end+1,1) = double(sym_Data_cell{i}); 
    
end   

end