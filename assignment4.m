%% Assignment 4

%% Question 1: Report on work done in PA 9

close all;
clear all;

% Parameters
R1 = 1;
Cap = 0.25;
R2 = 2;
L1 = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1; 
RO = 1000;

%setting up 8x8 matrix
G=zeros(8);
C=zeros(8);

%matrix C,G
G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 0 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0]; 

% Matrix G:
G

%Matrix C:
C

V1 = [];
V2 = [];

for Vin=-10:1:10
    Ffunc=[0; 0; 0; 0; 0; 0; 0; Vin];
    V=G\Ffunc;
    V1 = [V1 V(1)];
    V2 = [V2 V(5)];
end

figure(1)
hold on;
plot(-10:1:10, V1);
plot(-10:1:10, V2);
hold off;
legend('V3', 'VO');
title('DC Sweep');
ylabel('Voltage (V)');
xlabel('V (V)'); 

steps=1000;
Data=zeros(2,steps);
Data(1,:)=linspace(0,500,steps);
Vin=1;
for i=1:steps
    og=Data(1,i);
    Ffunc=[0; 0; 0; 0; 0; 0; 0; Vin];
    V=(G+1j*og*C)\Ffunc;
    Data(2,i)=V(5);
end

figure(2)
plot(Data(1,:),real(Data(2,:)));
title('AC plot - VO as a function of Omega');
ylabel('VO (V)');
xlabel('radians/s');

figure(3)
V2 = [];
stand = 0.05;
w = 3.14;
Vin = 1;
Ffunc=[0; 0; 0; 0; 0; 0; 0; Vin];

for w=1:1:10
    ep = (G+2*w^2*1j*C)\Ffunc;
    V2 = [V2 20*log10(abs(ep(5)/Ffunc(8)))];
end

semilogx(1:1:10, V2);
hold on;
title('AC Sweep');
xlabel('Radians/sec');
ylabel('dB');

figure(4);
cin = stand.*randn(5000,1) + Cap;
hold on;
V2 = [];
Ffunc=[0; 0; 0; 0; 0; 0; 0; Vin];

for index=1:5000
    C(1,1) = cin(index);
    C(2,1) = -cin(index);
    C(1,2) = -cin(index);
    C(2,2) = cin(index);
    ep = (G+2*pi*w*1j*C)\Ffunc;
    V2 = [V2 20*log10(abs(ep(5)/Ffunc(8)))];
end

title('Cap Sweep');
xlabel('Gain (dB)');
histogram(V2);

%% Question 2: Transient circuit simulation

iterations = 1000;
R1 = 1;
Cap = 0.25;
R2 =  2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1;
RO = 1000;
Cn = 0;
In = zeros(1,iterations);
time = 1;
dif = time/iterations;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0]; 

Vin = zeros(1,iterations);
Vin(0.03*iterations:iterations) = 1;
Ffunc = zeros(8,1,iterations);

for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

figure(5);
Vout = li(5,:,:);
Vout = Vout(1,:);
hold on;
plot(linspace(0,time,iterations), Vout);
plot(linspace(0,time,iterations), Vin);
title('Step input voltages');
legend('Vout', 'Vin');
xlabel('Time (s)');
ylabel('V (V)');

Vin = sin(linspace(0,1,iterations)*2*pi*1/0.03);
Ffunc = zeros(8,1,iterations);

for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);
figure(6);
hold on;
plot(linspace(0,time,iterations), Vout);
plot(linspace(0,time,iterations), Vin);
title('Voltage with sine input');
legend('Vout', 'Vin');
xlabel('Time (s)');
ylabel('Vo (V)');

Vin = gaussmf(linspace(0,1,iterations),[0.03 0.06]);

Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);
figure(7);
hold on;
plot(linspace(0,time,iterations), Vout);
plot(linspace(0,time,iterations), Vin);
hold off;
xlabel('Time (s)');
ylabel('Vout (V)');
title('Voltage for Gaussian Function');
legend('Vout', 'Vin');

Vin = sin(linspace(0,1,iterations)*2*pi*1/0.03);
Ffunc = zeros(8,1,iterations);

for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Ffunc = abs(fftshift(fft(Vout)));
figure(8);
hold on;
plot(((1:length(Ffunc))/iterations)-0.5,20*log10(Ffunc));

Ffunc = abs(fftshift(fft(Vin)));
plot(((1:length(Ffunc))/iterations)-0.5,20*log10(Ffunc));

xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Sine Function Frequency Response');

Vin = gaussmf(linspace(0,1,iterations),[0.03 0.06]);

Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Ffunc = abs(fftshift(fft(Vout)));

figure(9);
hold on;
plot(((1:length(Ffunc))/iterations)-0.5,20*log10(Ffunc));
Ffunc = abs(fftshift(fft(Vin)));
plot(((1:length(Ffunc))/iterations)-0.5,20*log10(Ffunc));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian Function Frequency Response');


%% Question 3: Circuit with noise

% New C matrix 

C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
C(7,:)=[0 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0]; 

C

Vin = gaussmf(linspace(0,1,iterations),[0.03 0.06]);
In = 0.001*rand(iterations,1);

Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);
figure(10);
hold on;
plot(linspace(0,time,iterations), Vout);
plot(linspace(0,time,iterations), Vin);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Voltages for Gassian plot with noise');
legend('Vout', 'Vin');

FF1 = abs(fftshift(fft(Vout)));
figure(11);
hold on;
plot(((1:length(FF1))/iterations)-0.5,20*log10(FF1));
FF1 = abs(fftshift(fft(Vin)));
plot(((1:length(FF1))/iterations)-0.5,20*log10(FF1));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
legend('Vout','Vin');
title('Gaussian pluse with noise');

figure(12);
hold on;
FF1 = abs(fftshift(fft(Vout)));
plot(((1:length(FF1))/iterations)-0.5,20*log10(FF1));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dBV)');
title('Voltages with noise with different Cn');

Cn = 0.1;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Cn = 0.01;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);
FF1 = abs(fftshift(fft(Vout)));
plot(((1:length(FF1))/iterations)-0.5,20*log10(FF1));

Cn = 0.0000000001;

G(1,:)=[1 -1 0 0 0 0 0 1];
C(1,:)=[Cap -Cap 0 0 0 0 0 0]; 
G(2,:)=[(-1/R1) (1/R2+1/R1) 0 0 0 1 0 0]; 
C(2,:)=[-Cap Cap 0 0 0 0 0 0];
G(3,:)=[0 0 1/R3 0 0 -1 0 0]; 
C(3,:)=[0 0 Cn 0 0 0 0 0]; 
G(4,:)=[0 0 0 alpha/R3 -1*alpha/R3 0 1 0]; 
C(4,:)=[0 0 0 0 0 0 0 0]; 
G(5,:)=[0 0 0 -1/R4 (1/R4+1/RO) 0 0 0]; 
C(5,:)=[0 0 0 0 0 0 0 0];
G(6,:)=[0 1 -1 0 0 0 0 0]; 
C(6,:)=[0 0 0 0 0 -L1 0 0]; 
G(7,:)=[0 0 -10 1 0 0 0 0];
C(7,:)=[0 0 0 0 0 0 0 0]; 
G(8,:)=[1 0 0 0 0 0 0 0]; 
C(8,:)=[0 0 0 0 0 0 0 0];

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);
FF1 = abs(fftshift(fft(Vout)));
plot(((1:length(FF1))/iterations)-0.5,20*log10(FF1));
legend('Cn = 0.1', 'Cn = 0.01', 'Cn = 0.0000000001');

%Changing Step Size

R1 = 1;
Cap = 0.25;
R2 =  2;
L = 0.2;
R3 = 10;
alpha = 100;
R4 = 0.1;
RO = 1000;
Cn = 0.00001;
iterations = 100;
time = 1;

dif = time/iterations;
Vin = gaussmf(linspace(0,1,iterations),[0.03 0.06]);
In = 0.001*rand(iterations,1);

Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);

figure(13);
hold on;
subplot(3,1,1);
plot(linspace(0,time,iterations), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 100');

iterations = 1000;
time = 1;
stand = 0.03;
mean = 0.06;
dif = time/iterations;
Vin = gaussmf(linspace(0,1,iterations),[stand mean]);
In = 0.001*rand(iterations,1);


Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end
Vout = li(5,:,:);
Vout = Vout(1,:);
subplot(3,1,2);
plot(linspace(0,time,iterations), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 1000');

iterations = 100000;
time = 1;
dif = time/iterations;
Vin = gaussmf(linspace(0,1,iterations),[stand mean]);
In = 0.001*rand(iterations,1);

Ffunc = zeros(8,1,iterations);
for i=1:iterations
    Ffunc(3,1,i) = -In(i);
    Ffunc(8,1,i) = Vin(i);
end

li = zeros(8,1, iterations);

for i=2:iterations
    ep = C/dif + G;
    li(:,:,i) = ep\(C*li(:,:,i-1)/dif + Ffunc(:,:,i));
end

Vout = li(5,:,:);
Vout = Vout(1,:);

subplot(3,1,3);
plot(linspace(0,time,iterations), Vout);
xlabel('Time (s)');
ylabel('Vout (V)');
title('Gaussian Function Transient Response with Various Time Steps - 100000');

%% Question 4: Non-linearity

%Another matrix would need to be used to represent the non-linear nature of
%the voltage.
