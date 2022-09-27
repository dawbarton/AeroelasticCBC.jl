% This Matlab Code computes Hopf point and Modal vectors of the grey-box
% model (using unsteady flutter model)
% from the free-decay response of the aerofoil
% Additionally it computes input form of Julia code that computes
% Centermanifold and Hopf normalform

% out1, out2,...,out6 gives nonlinear equation of motion in modal vectors
% r_out1, r_out2,...,r_out6 gives equation of motion in real vectors(heave, pitch, and so on)

clear all;
close all;
clc;

%% Importing the measured data (Free_decay)
load('Measured_data')
t_end=3.0; % Define the length of time data for system identification
% Description of variables (Stored in Measured_data.mat)
% data: measured response [heave;angle], b: wing semi chord, a: position of elastic axis to the semi-chord, m: mass ofthe wing
% m__T: total mass of the rig, w_a: modal freq of the pitch mode at U=0,
% w_h: modal freq of the heave mode at U=0, x__alpha: nondimensional
% distance between center of gravity and elastic axis, r_a: nondimensional radius of gyration 
% rho: air density, fs: sampling freuqency, 

% Free parameters: r_a, c__alpha,c__h
[r_a,c__alpha,c__h]=Lin_parameter_ID(data,b,a,m,m__T,w_a,w_h,x__alpha,r_a,c__alpha,c__h,Ts,fs,sen,t_end);

%%
I__alpha=m*(b*r_a)^2; % Mass moment of Inertia about the elastic axis
k__alpha=I__alpha*w_a^2; % pitch stifness
k__h=m__T*w_h^2; % heave stifness

U_max=30; % Maximum wind speed
wind_step=0.02; % Steps of wind speed to find Hopf point 
hopf_tolerance=1e-11; % zero-tolerance to find the zero real eigenvalue

[eig_value,modal_V,Uf]=Finding_flutter_speed(U_max,wind_step,hopf_tolerance,b,a,m,m__T,x__alpha,c__0,c__1,c__2,c__3,c__4,k__h,k__alpha,I__alpha,c__alpha,c__h,rho);


function [r_a,c__alpha,c__h] = Lin_parameter_ID(data,b,a,m,m__T,w_a,w_h,x__alpha,r_a,c__alpha,c__h,Ts,fs,sen,t_end)
h=data(2,:);
p=data(1,:);

h=h*sen; % Multiplying sensitivity of the measurement
zh=mean(h(end-100:end));  % Zero position of the displament (heave)
zp=mean(p(end-100:end)); % Zero position of the displament (angle)

time=1/fs:1/fs:1/fs*length(data(1,:)); 

[aa,cc]=max(abs(h));

sl=t_end*fs; % length of the time data

time_shift=1;
t_shift=fs*time_shift;

h1=h(cc+t_shift:cc+t_shift+sl)-zh;
p1=p(cc+t_shift:cc+t_shift+sl)-zp;

%% Parameter identification from the grey-box model (Using Matlab system identification toolbox)

in1=[h1;p1];
in1=transpose(in1);
data = iddata(in1,[],0.0002);

data.OutputName=[{'Heave (m)'}; {'Pitch (Rad)'}];
linear_model = idgrey('Flutter_rig', {x__alpha,r_a,c__alpha,c__h,w_a,w_h,m,m__T,b,a},'c');
linear_model.Structure.Parameters(1).Free = false;
linear_model.Structure.Parameters(5).Free = false;
linear_model.Structure.Parameters(6).Free = false;
linear_model.Structure.Parameters(7).Free = false;
linear_model.Structure.Parameters(8).Free = false;
linear_model.Structure.Parameters(9).Free = false;
linear_model.Structure.Parameters(10).Free = false;
opt = greyestOptions('InitialState','estimate','Display','off');
opt.Focus = 'stability';
linear_model = greyest(data,linear_model,opt);

r_a=linear_model.Structure.Parameters(2).Value; % Nondimensional Radius of gyration
c__alpha=linear_model.Structure.Parameters(3).Value;
c__h=linear_model.Structure.Parameters(4).Value;
end

function [J]=Flutter_Jac(U,b,a,m,m__T,x__alpha,c__0,c__1,c__2,c__3,c__4,k__h,k__alpha,I__alpha,c__alpha,c__h,rho)

MM = [b ^ 2 * pi * rho + m__T -a * b ^ 3 * pi * rho + b * m * x__alpha 0; -a * b ^ 3 * pi * rho + b * m * x__alpha I__alpha + pi * (0.1e1 / 0.8e1 + a ^ 2) * rho * b ^ 4 0; 0 0 1;];
DD = [c__h + 2 * pi * rho * b * U * (c__0 - c__1 - c__3) (1 + (c__0 - c__1 - c__3) * (1 - 2 * a)) * pi * rho * b ^ 2 * U 2 * pi * rho * U ^ 2 * b * (c__1 * c__2 + c__3 * c__4); -0.2e1 * pi * (a + 0.1e1 / 0.2e1) * rho * (b ^ 2) * (c__0 - c__1 - c__3) * U c__alpha + (0.1e1 / 0.2e1 - a) * (1 - (c__0 - c__1 - c__3) * (1 + 2 * a)) * pi * rho * (b ^ 3) * U -0.2e1 * pi * rho * (U ^ 2) * (b ^ 2) * (a + 0.1e1 / 0.2e1) * (c__1 * c__2 + c__3 * c__4); -1 / b a - 0.1e1 / 0.2e1 (c__2 + c__4) * U / b;];
KK = [k__h 2 * pi * rho * b * U ^ 2 * (c__0 - c__1 - c__3) 2 * pi * rho * U ^ 3 * c__2 * c__4 * (c__1 + c__3); 0 k__alpha - 0.2e1 * pi * (a + 0.1e1 / 0.2e1) * rho * (c__0 - c__1 - c__3) * (b ^ 2) * (U ^ 2) -0.2e1 * pi * rho * b * (U ^ 3) * (a + 0.1e1 / 0.2e1) * c__2 * c__4 * (c__1 + c__3); 0 -U / b c__2 * c__4 * U ^ 2 / b ^ 2;];

K1=-inv(MM)*KK;
D1=-inv(MM)*DD;

J1=[0,1,0,0,0,0];
J2=[K1(1,1),D1(1,1),K1(1,2),D1(1,2),K1(1,3),D1(1,3)];
J3=[0,0,0,1,0,0];
J4=[K1(2,1),D1(2,1),K1(2,2),D1(2,2),K1(2,3),D1(2,3)];
J5=[0,0,0,0,0,1];
J6=[K1(3,1),D1(3,1),K1(3,2),D1(3,2),K1(3,3),D1(3,3)];

J=[J1;J2;J3;J4;J5;J6];
end

function [eig_value,modal_V,Uf]=Finding_flutter_speed(U_max,wind_step,hopf_tolerance,b,a,m,m__T,x__alpha,c__0,c__1,c__2,c__3,c__4,k__h,k__alpha,I__alpha,c__alpha,c__h,rho)

U__n=wind_step:wind_step:U_max;
for ii=1:length(U__n)    
    U=U__n(ii);
    J=Flutter_Jac(U,b,a,m,m__T,x__alpha,c__0,c__1,c__2,c__3,c__4,k__h,k__alpha,I__alpha,c__alpha,c__h,rho);    
    value(ii,:)=eig(J);
    imaginary_part(ii,:)=imag(value(ii,:));
    [q1,q2]=sort(imaginary_part(ii,:)); 
    length_col=length(value(ii,:));    
    for jj=1:length_col
       value2(ii,jj)=value(ii,q2(jj));
    end    
    max_r_ev(ii)=max(real(value(ii,:)));
end

% Finding the location of wind speed where the sign of maximum real
% eigenvalue changes

for ii=1:length(U__n)-1    
    hopf_point(ii)= max_r_ev(ii)*max_r_ev(ii+1);
end

[aa,bb]=min(hopf_point);
max_r_ev2=min([abs(max_r_ev(bb)),abs(max_r_ev(bb+1))]);

% Finding the accurate Hopf Point which satisfies the tolerance
while abs(max_r_ev2)>hopf_tolerance
    U=(U__n(bb)+U__n(bb+1))/2;    
    %%
    J=Flutter_Jac(U,b,a,m,m__T,x__alpha,c__0,c__1,c__2,c__3,c__4,k__h,k__alpha,I__alpha,c__alpha,c__h,rho);
    
    [eig_value]=eig(J);
    [modal_V,Diagonal]=eig(J);   
    max_r_ev2=max(real(eig_value));    
    if max_r_ev2>0
        U__n(bb+1)=U;
    else
        U__n(bb)=U;
    end
end
Uf=U;
end



