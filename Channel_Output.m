%Most of the conventions used are similar to those in "Nested Sparse Approximation: Structured Estimation of V2V Channels Using Geometry-Based Stochastic Channel Model"
K = 50; 
M = 180; 
Nr = 100; 
%reading the data from "h_t.txt" which has details about the channel impulse response
h_n_m = dlmread('h_t-v5.txt');
H_k_m = zeros(2*K+1,M);
%finding DFT
for k=-K:K
    for m=0:M-1
        for n=0:Nr-1
            H_k_m(k+K+1,m+1)= H_k_m(k+K+1,m+1)+ (1/(2*K+1))*(h_n_m(n+1,m+1)*(cos(2*pi*n*k/(2*K+1))-sin(2*pi*n*k/(2*K+1))*sqrt(-1)));
        end
    end
end
s_n = 1.5*(randi([1 2],1,Nr+M-1)*2 - 3); % transmitted pilot sequence %dlmread('ip-v5.txt'); 
N = (2*K+1)* M;
x = zeros(N,1);
j=1;
for m=1:M
    for k=1:2*K+1
        x(j,1)= H_k_m(k,m);
        j=j+1;
    end
end

S = zeros(Nr,N); %input block data matrix of size Nr x N

 % Creating a vandermonde matrix of size Nr x (2K+1)
omega = zeros(Nr, 2*K+1);
S_m = zeros(Nr, Nr); %size Nr x Nr

%Updating the vandermonde matrix
for i=0:Nr-1
    for j=0:2*K
        omega(i+1,j+1)=(cos(2*pi/(2*K+1))+sin(2*pi/(2*K+1))*sqrt(-1))^(i*(j-K));
    end
end
p=1;

%defining S matrix
for m=0:M-1
    for n=0:Nr-1
        S_m(n+1,n+1) = s_n(-m+n+M);
    end
    S_m2 = S_m * omega;
    %Copying S_m2 into S
    for o=1:Nr
        for q=1:Nr
            S(q,p)=S_m2(q,o);
        end
        p=p+1;
    end
end
y_n = S*x;
y_n_mod = real(y_n) ;%abs(y_n);
figure(1);
bar(y_n_mod);
title('Output samples');
ylabel('y[n]]');
figure(2);
bar(s_n);
h_t_mod = abs(h_n_m);
title('Input samples');
xlabel('t (x10^-9 s)'), ylabel('s[n]');
figure(3);
bar3(h_t_mod,0.01);
title('Discrete time impulse response');
ylabel('t (x10^-9 s)'), xlabel('\tau (0.01x10^-6s)'), zlabel('h[n,m]');
x_mod = abs(x);
figure(4);
bar(x_mod,0.01);

%writing the data into a file
dlmwrite('S-v4.txt',S);
