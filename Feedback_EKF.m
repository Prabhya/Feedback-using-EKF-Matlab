%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating Feedback control using ekf.
%
% - Prabhjeet Singh Arora
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

% Time
dt = 0.1;
tf = 60;
t  = 0:dt:tf;
m  = length(t);

% System m = 1, c = 2
% x_dot = [x_2;2x_1 -3x_2]+ BKx + w_k
A = [0 1;-2 2];
B = [0;1];
k1 = 3;
k2 = 5;
FB = [k1 k2];
C = [1 0];
G = [0;1];

x0 = [5;-10];
x_true = zeros(2,m);
x_true(:,1)=x0;

xact_store = zeros(2,m);
yact_store = zeros(1,m);
xe_store = zeros(2,m);

Q = 10^(-1);
R = 10^(-2);

Qt = sqrt(Q)*randn(1,m);
x0act = x0;
y0act = C*x0act + sqrt(R)*randn;
yact_store(1) = y0act;
xact_store(:,1) = x0act; 

pcov0 = eye(2);
pcov = pcov0;
pcov_store = zeros(2,m);


K = pcov*C'/(C*pcov*C'+R);
xe0 = x0;
xe = xe0;
ye = C*xe;
%xe = xe+K*(yact_store(1)-ye);
%pcov = (eye(2)-K*C)*pcov;
xe_store = zeros(2,m);
xe_store(:,1) = xe;
pcov_store(:,1) = diag(pcov);
q= Q;

%Qe = sqrt(Q)*randn(1,m);
for i =2:m
    Xt = ode4(@fun,[t(i-1) t(i)],[x_true(:,i-1);x_true(:,i-1);0],A,B,FB,G);
    x_true(:,i) = Xt(end,1:2)';
    
    X = ode4(@fun,[t(i-1) t(i)],[xact_store(:,i-1);xe;Qt(i-1)],A,B,FB,G);
    xact_store(:,i) = X(end,1:2)';
    yact_store(i) = C*xact_store(:,i)+ sqrt(R)*randn;

    X1 = ode4(@fun2,[t(i-1) t(i)],[xe;xe;0;reshape(pcov,[],1)],A,B,FB,G,q);
    xe = X1(end,1:2)';
    pcov = reshape(X1(end,6:end)',2,2);
    ye = C*xe;
    K = pcov*C'/(C*pcov*C'+R);
    xe = xe+K*(yact_store(i)-ye);
    pcov = (eye(2)-K*C)*pcov;
    xe_store(:,i) = xe;
    pcov_store(:,i) = diag(pcov);
    %if norm(xe)>200
    %    disp(i)
    %    break
    %end
end
%sig = 3*pcov_store.^0.5;

figure(1)
plot(t,xact_store(1,:),'.-r',t,x_true(1,:),'.b')
%axis([0 60 -0.5 0.5])
grid on
title('Creating state-feedback system using EKF - Position',"FontSize",18)
xlabel('Time (sec)',"FontSize",18)
ylabel('Position (m)',"FontSize",18)
legend('Actual state-feedback using EKF','True feedback system w/o noise')

figure(2)
plot(t,xact_store(2,:),'.-r',t,x_true(2,:),'.b')
grid on
title('Creating state-feedback system using EKF - Velocity',"FontSize",18)
xlabel('Time (sec)',"FontSize",18)
ylabel('Velocity (m s^-^1)',"FontSize",18)
legend('Actual state-feedback using EKF','True feedback system w/o noise')



function f = fun(~,x,A,B,FB,G)
    x_state = x(1:2);
    x_fb = x(3:4);
    q = x(5);
    x_dot = A*x_state -B*FB*x_fb+G*q;
    f = [x_dot;x_fb;q];
end

function f = fun2(~,x,A,B,FB,G,Q)
    x_state = x(1:2);
    x_fb = x(3:4);
    q = x(5);
    pcov = reshape(x(6:end),2,2);
    x_dot = A*x_state-B*FB*x_fb;
    pcov_dot = A*pcov+pcov*A'+G*Q*G';
    f=[x_dot;x_fb;q;reshape(pcov_dot,[],1)];
end