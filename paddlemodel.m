%% Paddle Model
mah
function paddlemodel()
global m M g l

% prep fig
clear
clc
figure(1);clf;
subplot 221;hold on
subplot 223;hold on
subplot 122;hold on

subplot 122; plot(0,0,'+', 'DisplayName', 'origin')

% parameters
m = 0.2;
M = 1.0;
g = 9.8;
l = 1.0;

params = [m M g l];

BK=[0 0];
tspan=[0 30];

% system definition
A = [
    0,0,1,0;
    0,0,0,1;
    0,-g,0,0;
    -12*g*m/(l^2*M),0,0,0];
B = [0;0;0;12/(l^2*M)];
C = [1,0,0,0;0,1,0,0];
D = 0;

plant = ss(A,B,C,D);

% sim definition
t=0:.1:6;
tspan=[0 30];
x0=[
    0.0 0.0 0.0 0.0;
    0.1 0.0 0.0 0.0;
    0.0 -deg2rad(10.0) 0.0 0.0;
    -0.2 deg2rad(20.0) 0.0 0.0;
    0.3 deg2rad(30.0) 0.0 0.0;
    ]';
x0_hat = x0; %start out with good estimate

n_test = size(x0, 2);

% for i=1:2,
%     [~,t,x]=lsim(plant,zeros(size(t)),t,x0(:,i));
%     subplot 221;plot(t,x(:,1),'g--','DisplayName','lin plant')
%     subplot 223;plot(t,x(:,2),'g--','DisplayName','lin plant')
%     subplot 122;plot(x(:,1),x(:,2),'g--','DisplayName','lin plant')
% end

%% Controller
%full state, nonlinear
Q = diag([1.0,1.0,0.0,0.0]);
R = 0.2;


% lqr ...
[Kfull,~,~] = lqr(A,B,Q,R);
Kfull = -Kfull;

% pole placement ...
s1 = -1 / 0.05;
s2 = -1 / 0.05;
s3 = -1 / 1.0;
s4 = -1 / 1.0;

k1 = m*g + l^2*M*s1*s2*s3*s4/(12*g);
k2 = -1/12*l^2*M*(s3*s4+s2*(s3+s4)+s1*(s2+s3+s4));
k3 = -l^2*M*(s1*s2*s3+s2*s3*s4+s1*(s2+s3)*s4)/(12*g);
k4 = 1/12*l^2*M*(s1+s2+s3+s4);

Kpole = [k1 k2 k3 k4];

Kpole
Kfull

BK=B*Kpole;

% observer ...

% Plots

% lin, cl
sys_cl=ss(A+BK,B,C,D);
for i=1:n_test,
    [~,t,x]=lsim(sys_cl,zeros(size(t)),t,x0(:,i));
    subplot 221;plot(t,x(:,1),'r--','DisplayName','lin full')
    subplot 223;plot(t,x(:,2),'r--','DisplayName','lin full')
    subplot 122;plot(x(:,1),x(:,2),'r--','DisplayName','lin full')
end

% non-lin, cl
for i=1:n_test,
    [t,x] = ode45(@(t,x) odefun(t,x,BK, params), tspan, x0(:,i));
    subplot 221;plot(t,x(:,1),'r','DisplayName','full')
    subplot 223;plot(t,x(:,2),'r','DisplayName','full')
    subplot 122;plot(x(:,1),x(:,2),'r','DisplayName','full')
end

% non-lin, cl + non-lin obs
% L=[-3 -4]';
% x0t=[x0(:,1)' 1*x0(:,1)'-[1 0];x0(:,2)' 1*x0(:,2)'-[1 0]]';
% for i=1:n_test,
%     [t,x] = ode45(@(t,x) odefun2(t,x,BK,L*C,xeq), tspan, x0t(:,i));
%     subplot 221;plot(t,x(:,1)-x(:,3),'c','DisplayName','estimate','linewidth',2)
%     subplot 223;plot(t,x(:,2),'c','DisplayName','estimate','linewidth',2)
%     subplot 122;plot(x(:,1)-x(:,3),x(:,2)-x(:,4),'c','DisplayName','estimate','linewidth',2)
% 
%     subplot 221;plot(t,x(:,1),'b','DisplayName','observer')
%     subplot 223;plot(t,x(:,2),'b','DisplayName','observer')
%     subplot 122;plot(x(:,1),x(:,2),'b','DisplayName','observer')
% end

grid;legend('show')
end

%% System State Transition
function fx = f(s, p)
tmp = num2cell(p);
[m, M, g, l] = deal(tmp{:});

x = s(1);
t = s(2);
v = s(3);
w = s(4);

ax = -g * sin(t);
at = -g * m * x * cos(t) / ...
    (l^2*M/12+m*x^2);

fx = [v,w,ax,at]';
end

%% ODE v1
function dxdt = odefun(~,x,BK,p)
dxdt = f(x,p)+BK*x;
end