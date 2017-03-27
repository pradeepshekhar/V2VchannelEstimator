K = 50; 
M = 180;
Nr = 100;
N = (2*K+1)* M;
y_n = dlmread('op-v5.txt');
S = dlmread('S-v5.txt');
cvx_begin 
    variable x(N)
    minimize( norm( x, 1) )
    subject to
        norm((y_n - S*x),2) <= 0.005;%abs(y_n - S*x) <= 0.005;
cvx_end
%x_mod = abs(x) ;
%figure(1);
%bar(x_mod);
