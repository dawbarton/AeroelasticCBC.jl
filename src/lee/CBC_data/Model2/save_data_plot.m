clear all;
close all;

% Description of variables
% amp: amplitude of the noninvasive control target measured in the experiment
% U: Wind speed, freq: frequency of the oscillation, error1: (control error)^2
% error_t: control error computed from the time-series in data

load('CBC_v20_7') 
sen=18/100; % sensitivity of the heave response

amp(1)=amplitude; % amplitude of the heave response
U(1)=Wind_sped;
freq(1)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(1)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(1)=ah;
error_t(1)=at-ah;
%%
load('CBC_v21_6')%%

amp(2)=amplitude; % amplitude of the heave response
U(2)=Wind_sped;
freq(2)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(2)=ah;
error_t(2)=at-ah;
%%

load('CBC_v22_4')

amp(3)=amplitude; % amplitude of the heave response
U(3)=Wind_sped;
freq(3)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(3)=ah;
error_t(3)=at-ah;

load('CBC_v23_3') %%

amp(4)=amplitude; % amplitude of the heave response
U(4)=Wind_sped;
freq(4)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(4)=ah;
error_t(4)=at-ah;

load('CBC_v24_1')
amp(5)=amplitude; % amplitude of the heave response
U(5)=Wind_sped;
freq(5)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(5)=ah;
error_t(5)=at-ah;

load('CBC_v25_0')
amp(6)=amplitude; % amplitude of the heave response
U(6)=Wind_sped;
freq(6)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(6)=ah;
error_t(6)=at-ah;


load('CBC_v25_8')
amp(7)=amplitude; % amplitude of the heave response
U(7)=Wind_sped;
freq(7)=freq_actual;

ma=mean(data(4,:)); % mean position of the angle
a=data(4,:)-ma;
angle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);

mt=mean(data(2,:));
tt=data(2,:)-mt;
at=max(tt);

Uamp2(7)=ah;
error_t(7)=at-ah;
%% stable data

load('CBC_stable_v20_7')
s_amp(1)=amplitude;
s_U(1)=Wind_sped;
s_freq(1)=freq_actual;

ma=mean(data(4,:));
a=data(4,:)-ma;
sangle(1)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);
Samp2(1)=ah;

load('CBC_stable_v21_5')
s_amp(2)=amplitude;
s_U(2)=Wind_sped;
s_freq(2)=freq_actual;

ma=mean(data(4,:));
a=data(4,:)-ma;
sangle(2)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);
Samp2(2)=ah;

load('CBC_stable_v22_4')
s_amp(3)=amplitude;
s_U(3)=Wind_sped;
s_freq(3)=freq_actual;

ma=mean(data(4,:));
a=data(4,:)-ma;
sangle(3)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);
Samp2(3)=ah;

load('CBC_stable_v23_2')
s_amp(4)=amplitude;
s_U(4)=Wind_sped;
s_freq(4)=freq_actual;

%%
ma=mean(data(4,:));
a=data(4,:)-ma;
sangle(4)=max(a);

mh=mean(data(1,:));
h=data(1,:)-mh;
ah=max(h);
Samp2(4)=ah;

%% Flutter speed and frequency estimation

amp=amp*sen;
s_amp=s_amp*sen;
Uamp2=Uamp2*sen;
Samp2=Samp2*sen;

amp2=amp.^2;

p = polyfit(U,amp2,1);

uu=20:0.01:27;

reg_curve=p(1)*uu+p(2);

figure(1)
stem(U,amp2,'LineStyle','none');
hold on;
plot(uu,reg_curve);
xlabel('Wind Speed(m/sec)')
ylabel('Square of heave amplitude (m^2)')
title('Flutter speed estimation')

flutter_speed=-p(2)/p(1);

delta=flutter_speed-U;

p2 = polyfit(delta,freq,1);

delta2=0:0.01:7;
reg_curve_f=p2(1)*delta2+p2(2);

flutter_freq=p2(2);

figure(2)
stem(delta,freq,'LineStyle','none')
hold on

plot(delta2,reg_curve_f)
xlabel('delta (m/sec)')
ylabel('Frequency (Hz)')
ylim([2.36 2.42])
title('Flutter Frequency estimation')


figure(3)
stem(U,amp,'LineStyle','none');
hold on;
stem(s_U,s_amp,'LineStyle','none');
xlabel('Wind Speed (m/sec)')
ylabel('Heave Amplitude (m)')
legend('Unstable LCO','Stable LCO')
title('Bifurcation diagram plotted by noninvasive control target')

figure(4)
stem(U,freq,'LineStyle','none');
hold on;
stem(s_U,s_freq,'LineStyle','none');
ylim([2.38 2.45])
xlabel('Wind Speed (m/sec)')
ylabel('Frequency(Hz)')
legend('Unstable LCO','Stable LCO')
title('Frequency of LCOs')

figure(5)
stem(U,Uamp2,'LineStyle','none');
hold on;
stem(s_U,Samp2,'LineStyle','none');
xlabel('Wind Speed (m/sec)')
ylabel('Heave Amplitude (m)')
legend('Unstable LCO','Stable LCO')
title('Bifurcation diagram plotted by time-series')
