%Most of the conventions used are similar to those in "A Geometry-Based Stochastic MIMO Model for Vehicle-to-Vehicle Communications"
%Total time frame is divided into 
t_res = 10^-9;
tau_quant = 0.01*10^-6;
t_samples = 100; %Nr - 200
tau_samples = 180; %M - 256
Co = 3*10^8; %propagation speed in m/s
lambda_v = 0.0508 ; %considering frequency of opereation = 5.9 GHz
x_min = 0;
x_max = 500; % 1000;
w = 18; %width of road
y1_di = 6.5;
y2_di = -6.5;
y1_sd = 6.5;
y2_sd = -6.5;
var_sd = 4;
W_di = 5;
N_total = 430;
N_md = 10;
N_sd = 20;
N_di = 400;
N_lanes = 4;
v_p = 19.444; %velocity of md scatterers : 70kmph
v_r = 22.222; %velocity of Rx : 80kmph
v_t = 27.78; %velocity of Tx : 100kmph
x_t = 200; %10 - intial position of Tx
y_t =-1;
x_r = 300; %200 - intial position of Rx
y_r =1;
var_glos = 1; %variance corresponding to stochastic amplitude gain of scatterers
var_gsd = 1;        
var_gmd = 1;
x_md = x_min + (x_max - x_min)*rand(1,N_md);
x_sd = x_min + (x_max - x_min)*rand(1,N_sd);
x_di = x_min + (x_max - x_min)*rand(1,N_di);
y_sd1 = (sqrt(var_sd)) * randn(1,N_sd/2) + y1_sd;
y_sd2 = (sqrt(var_sd)) * randn(1,N_sd/2) + y2_sd;
y_sd = [y_sd1 y_sd2];
y_di1 = (y1_di - W_di/2) + W_di * rand(1,N_di/2);
y_di2 = (y2_di - W_di/2) + W_di * rand(1,N_di/2);
y_di = [y_di1 y_di2];
y_md = randi([1 N_lanes],1,N_md)*2 - 5;
tau = zeros(1,N_total+1);
nu = zeros(1,N_total+1);
h_t = zeros(t_samples,tau_samples);
k = 1;
d_ref = 1; % 2;
%LOS parameters
Go_los = 0.562341; %-5dB
n_los = 1.8;
g_los = (sqrt(var_glos))*randn(1,1); %stochastic amplitude gain

phi_sd = 2*pi*rand(1,N_sd) ;%phase shift of static discrete scatterers
phi_md = 2*pi*rand(1,N_md);%phase shift of moving discrete scatterers
n_sd = 3.5*rand(1,N_sd); 
n_md = 3.5*rand(1,N_md); 

%************************************************
Go_sd = -89+24*n_sd;%in dB.... to be changed?
Go_md = -89+24*n_md; %in dB.... to be changed?
for n=1:N_sd
    Go_sd(n) = 10^(Go_sd(n)/20);
end
for n=1:N_md
    Go_md(n) = 10^(Go_md(n)/20);
end
%************************************************

g_sd = (sqrt(var_gsd))*randn(1,N_sd); %stochastic amplitude gain
g_md = (sqrt(var_gmd))*randn(1,N_md); %stochastic amplitude gain
 
var_amp_di = 1;%variance of distribution representing the diffuse scatterer amplitudes **************
Go_di = 14.1253; %23dB
n_di = 5.4; %pathloss exponent - value taken from paper mentioned in first line
Cr = (sqrt(var_amp_di)) * randn(1,N_di); %distribution representing the diffuse scatterer amplitudes
t=0;
E= 0.0001 ;%parameter used to increment the slowly varying stochastic amplitude gains
%varying part of the stochastic amplitude gain
v_los = randn(1,1); 
v_sd = randn(1,N_sd);
v_md = randn(1,N_md);

for l=1:t_samples   
    %evaluating contribution from LOS component
    d0 = sqrt((y_r-y_t)*(y_r-y_t)+(x_r-x_t)*(x_r-x_t));
    tau(k) = d0/Co ;
    tau_b = tau(k)/tau_quant;
    tau_c = fix(tau_b);
    if (tau_b - tau_c)<= 0.5
        tau_f = tau_c;
    else
        tau_f = tau_c+1;
    end  
    nu(k) = ((v_t-v_r)*((x_r-x_t)/d0))/lambda_v;
    amp = (g_los)*(Go_los^0.5)*((d_ref/d0)^(n_los/2))*(cos(2*pi*nu(k)*t)-sin(2*pi*nu(k)*t)*sqrt(-1));
    h_t(l,tau_f)=h_t(l,tau_f)+amp;
    k=k+1;
    %evaluating contributions from Mobile discrete scatterers
    for i=1:N_md
        d0=sqrt((y_md(i)-y_t)*(y_md(i)-y_t)+(x_md(i)-x_t)*(x_md(i)-x_t));
        d1=sqrt((y_md(i)-y_r)*(y_md(i)-y_r)+(x_md(i)-x_r)*(x_md(i)-x_r));
        tau(k) = (d0 + d1)/Co;
        tau_b = tau(k)/tau_quant;
        tau_c = fix(tau_b);
        if (tau_b - tau_c)<= 0.5
            tau_f = tau_c;
        else
            tau_f = tau_c+1;
        end        
        nu(k) = (((v_t-v_p)*((x_md(i)-x_t)/d0))+((v_r-v_p)*((x_md(i)-x_r)/d1)))/lambda_v;
        amp = (g_md(i))*(cos(phi_md(i))-sin(phi_md(i))*sqrt(-1))*(Go_md(i)^0.5)*((d_ref/(d0+d1))^(n_md(i)/2))*(cos(2*pi*nu(k)*t)-sin(2*pi*nu(k)*t)*sqrt(-1));
        h_t(l,tau_f) = h_t(l,tau_f)+amp;        
        k=k+1;
    end
    %evaluating contributions from static discrete scatterers
    for i=1:N_sd
        d0 =sqrt((y_sd(i)-y_t)*(y_sd(i)-y_t)+(x_sd(i)-x_t)*(x_sd(i)-x_t));
        d1 = sqrt((y_sd(i)-y_r)*(y_sd(i)-y_r)+(x_sd(i)-x_r)*(x_sd(i)-x_r));
        tau(k) = (d0 + d1)/Co;
        tau_b = tau(k)/tau_quant;
        tau_c = fix(tau_b);
        if (tau_b - tau_c)<= 0.5
            tau_f = tau_c;
        else
            tau_f = tau_c+1;
        end        
        nu(k) = (((v_t)*((x_sd(i)-x_t)/d0))+((v_r)*((x_sd(i)-x_r)/d1)))/lambda_v;
        amp = (g_sd(i))*(cos(phi_sd(i))-sin(phi_sd(i))*sqrt(-1))*(Go_sd(i)^0.5)*((d_ref/(d0+d1))^(n_sd(i)/2))*(cos(2*pi*nu(k)*t)-sin(2*pi*nu(k)*t)*sqrt(-1));
        h_t(l,tau_f)=h_t(l,tau_f)+amp;                
        k=k+1;
    end
    %evaluating contributions from diffuse scatterers
    for i=1:N_di
        d0 = sqrt((y_di(i)-y_t)*(y_di(i)-y_t)+(x_di(i)-x_t)*(x_di(i)-x_t));
        d1 = sqrt((y_di(i)-y_r)*(y_di(i)-y_r)+(x_di(i)-x_r)*(x_di(i)-x_r));
        tau(k) = (d0 + d1)/Co;
        tau_b = tau(k)/tau_quant;
        tau_c = fix(tau_b);
        if (tau_b - tau_c)<= 0.5
            tau_f = tau_c;
        else
            tau_f = tau_c+1;
        end        
        nu(k) = (((v_t)*((x_di(i)-x_t)/d0))+((v_r)*((x_di(i)-x_r)/d1)))/lambda_v;
        amp = ((Go_di^0.5)*(Cr(i))*((d_ref/(d0*d1))^(n_di/2)))*(cos(2*pi*nu(k)*t)-sin(2*pi*nu(k)*t)*sqrt(-1));
        h_t(l,tau_f)=h_t(l,tau_f)+amp;                
        k=k+1;
    end

    x_t = x_t + v_t*t_res;
    x_r = x_r + v_r*t_res;
    for m=1:N_md
        x_md(m)=x_md(m)+ v_p*t_res;
    end
    k=1;
    t=t+t_res;
    %incrementing the slowly varying stochastic amplitude gains
    g_los = (sqrt(1-E))* g_los + sqrt(E) * v_los ;
    for i=1:N_sd
        g_sd(i) = (sqrt(1-E))* g_sd(i) + sqrt(E) * v_sd(i)  ;
    end
    for i=1:N_md
        g_md(i) = (sqrt(1-E))* g_md(i) + sqrt(E) * v_md(i) ;
    end
end
h_t_mod = abs(h_t);
figure(1);
bar3(h_t_mod,0.01); %mesh, surf
title('Discrete time impulse response');
xlabel('t (x10^-9 s)'), ylabel('\tau (0.01x10^-6s)'), zlabel('h[n,m]');
figure(2);
plot(tau,nu,'*');
title('Delay Doppler representation of the V2V channel');
xlabel('\tau (s)'), ylabel('\nu (Hz)');
%h_t_phi = angle(h_t);
%figure(3);
%bar3(h_t_phi,0.1);
%dlmwrite('h_t-v2.txt',h_t)