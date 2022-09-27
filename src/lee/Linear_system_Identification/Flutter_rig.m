function [A,B,C,D] = Flutter_rig(x__alpha,r_a,c__alpha,c__h,w_a,w_h,m,m__T,b,a,Ts)
% Equation of motion of Flutter equation
% Input of system identification toolbox

rho=1.204; % Air density
I__alpha=m*(b*r_a)^2; % Mass moment of Inertia

k__alpha=I__alpha*w_a^2;
k__h=m__T*w_h^2;

MM = [m__T+pi*rho*b^2,m*x__alpha*b-a*pi*rho*b^3;m*x__alpha*b-a*pi*rho*b^3,I__alpha+pi*(1/8+a^2)*rho*b^4];
DD = [abs(c__h),0;0,abs(c__alpha)];
KK = [k__h,0;0,k__alpha];

II=eye(2);
OO=zeros(2,2);

A=[OO,II;-inv(MM)*KK,-inv(MM)*DD];
B=zeros(4,0);
C=[1,0,0,0;0,1,0,0];
D=zeros(2,0);
end
