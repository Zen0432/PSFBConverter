clc;clear;close all;

R_i = 100;
R_ct = 0;
C_dl = 1;
V_oc = 0;

n = 1/4;
L_lk = 5e-6;
f_s = 20000;
V_d = 400;
L = 100e-6;
C = 100e-6;
R_eq = 4*L_lk*f_s*n^2;

d = 0.5;

A = [ -R_eq/L    -1/L      ;
       1/C       -1/(R_i*C);]      
B=[2*n*V_d/L ; 0 ];
C=[0 1];

Co=ctrb(A,B)
rank(Co)
ob=obsv(A,C)
Ro=rank(ob)

t=0:0.00001:0.04;

u=d*ones(size(t));

sys_org = ss(A,B,C,0)
eig(sys_org)
tf(sys_org)


[y,t] = lsim(sys_org,u,t);
figure(1)
plot(t,y)
hold on
title("Original System Response Duty Cycle 0.5(Max duty cycle)")
xlabel("Time")
ylabel("Amplitude")
hold off;

K=place(A,B,11.6*[(-500+700i),(-500-700i)])

sys_mod = ss(A-B*K,B,C,0)
eig(sys_mod)
tf(sys_mod)

t1=0:0.00001:0.005;
u1 = d*ones(size(t1));

[y,t,x] = lsim(sys_mod,u1,t1);

figure(2)
plot(t,y)
hold on
title("System Response After Addition of Gains Matrix K Duty Cycle 0.5(Max duty cycle)")
xlabel("Time")
ylabel("Amplitude")
hold off;


figure(3)
hold on;
bode(sys_org)
bode(sys_mod)
hold off;

L=place(A',C',116*[(-500+700i),(-500-700i)])'


At=[A-B*K B*K;zeros(size(A)) A-L*C];
Bt=[B; zeros(size(B))];
Ct=[C zeros(size(C))];
sys_o=ss(At,Bt,Ct,0);

%figure(5)
% [y,t] = lsim(sys_o,u1,t1,[[0,0] [0,0]]);
% plot(t,y)
% hold on
% title("Observer Response")
% xlabel("Time")
% ylabel("Amplitude")
% hold off;

n=2;
[y,t1,x]=lsim(sys_o,u1,t1,[[0,0] [0,0]]);
e=x(:,n+1:end);
x=x(:,1:n);
xe=x-e;

x1=x(:,1);x2=x(:,2);
xe1=xe(:,1);xe2=xe(:,2);
figure(4)
plot(t1,x1,'-r',t1,xe1,'--r',t1,x2,'-b',t1,xe2,'--b')
hold on
title("Observer Response(Blue) and Error");
xlabel("Time");
ylabel("Amplitude");
hold off;


