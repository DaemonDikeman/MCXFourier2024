%Power Spectrum
clf
clear

%frequency resolution desired: 0.11
%FR = fs/fft size
%0.10 = 7/fft -> fft size = 7/.10 = ~70
%current resolution: 0.0025 = 7/x -> fft size of 2800


v=36;
m=.01;
fr=0.11;

fs=700*m;
T=1/fs;
t=0:T:(T*2*round(fs/fr))-T;

x=linspace(10*m,66*m, v);
z=linspace(10*m,66*m, v);
adjust=zeros(1,v);

for i=1:v
    adjust(i)=0.2*x(i)*(-0.5+rand);
    z(i)=z(i)+adjust(i);
end
yn=cell(v,1);
yn1=cell(v,1);
sn = [1.3,1.38,1.46,1.54,1.62,1.7,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.10...
    1.06,1.02,0.98,0.94,0.90,0.86,0.82,0.76,0.7,0.64,0.58,0.52,0.46,0.4...
    0.375,0.35,0.325,0.3,0.275,0.25,0.225,0.2];
sn = sn.^(2);

for i=1:v
   yn{i}=sn(i)*sin(z(i)*pi*(t+0.50/z(i) - 0.5));
   yn1{i}=sn(i)*sin(x(i)*pi*(t+0.50/z(i) - 0.5));
end

yavg=zeros(1,length(t));
yavg1=zeros(1,length(t));

for i=1:length(t)
    y=0;
    y1=0;
    for j=1:length(yn)
        y=y+yn{j,1}(i);
        y1=y1+yn1{j,1}(i);
    end
    yavg(i)=y/v;
    yavg1(i)=y1/v;
end
L=length(yn{1,1});

yx=fft(yavg)./numel(yavg);
yx1=fft(yavg1)./numel(yavg1);
xpow=abs(yx(1:L/2+1));
xpow1=abs(yx1(1:L/2+1));

fn=cell(v,1);
fn1=cell(v,1);
%n=2^nextpow2(L);
for i=1:v
    ym=fft(yn{i,1})./numel(yn{i,1});
    ym1=fft(yn1{i,1})./numel(yn{i,1});
    power=abs(ym(1:L/2+1));
    power1=abs(ym1(1:L/2+1));
    %power=abs(ym(1:L/2)).^2;
    %P2=abs(ym);
    %P1=P2(1:L/2+1);
    %P1(2:end-1)=2*P1(2:end-1);
    fn{i,1}=power;
    fn1{i,1}=power1;
%     ym(1) = [];
     
%     fn{i,1}=power;
%     fn{i,1}(1)=[];
end
nyquist=1;
freq=fs*(0:L/2)/(L/2)*nyquist;
%period=1./freq;
%freq=2*(1:(L/2))/L;

pavg=zeros(1,length(power));
pavg1=zeros(1,length(power1));
for i=1:length(pavg)-1
    p=0;
    p1=0;
    for j=1:v
        p=p+fn{j,1}(i);
        p1=p1+fn1{j,1}(i);
    end
    pavg(i)=p/v;
    pavg1(i)=p1/v;
end

%d=pspectrum(yn{15});


figure(1)
subplot(2,1,1)
plot(t, yn{1})
hold on
for i=2:v
    plot(t, yn{i});
end
hold off
grid
xlim([0,10])
subplot(2,1,2)
plot(t, yn1{1})
hold on
for i=2:v
    plot(t, yn1{i});
end
hold off
grid
xlim([0,10])
title('All Sinusoids')

figure(2)
subplot(2,1,1)
plot(t, yavg)

grid
subplot(2,1,2)
plot(t, yavg1)
grid
title('Average signal')

figure(3)
subplot(2,1,1)
plot(freq, 6*pavg, '-+')
xlim([8*m,80*m])
ylim([25*m, 200*m])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8])
yticks([0 0.4 0.8 1.2 1.6 2.0])
grid
title('Power Spectrum')
subplot(2,1,2)
plot(freq, 6*pavg1, '-+')
xlim([8*m,80*m])
ylim([25*m, 200*m])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8])
yticks([0 0.4 0.8 1.2 1.6 2.0])
grid
title('Power Spectrum')

figure(4)
subplot(2,1,1)
plot(freq, 10*xpow, '-+')
xlim([10*m,80*m])
ylim([25*m, 200*m])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8])
yticks([0 0.4 0.8 1.2 1.6 2.0])
grid
title('Power Spectrum')
subplot(2,1,2)
plot(freq, 10*xpow1, '-+')
xlim([10*m,80*m])
ylim([25*m, 200*m])
xticks([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8])
yticks([0 0.4 0.8 1.2 1.6 2.0])
grid
title('Power Spectrum')

% figure(5)
% subplot(2,1,1)
% plot(log(freq), log(pavg))
% grid
% subplot(2,1,2)
% plot(log(freq), log(pavg1))
% grid