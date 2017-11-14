function nonlinobserver()
clear
clc

figure(1);clf;
subplot 221;hold on
subplot 223;hold on
subplot 122;hold on

BK=[0 0];
xeq=[0 0]';

tspan=[0 30];
x0=[-1 1.327;3 -1.326]';
for i=1:2,
     [t,x] = ode45(@(t,x) odefun(t,x,BK,xeq), tspan, x0(:,i));
     subplot 221;plot(t,x(:,1),'g','DisplayName','plant')
     subplot 223;plot(t,x(:,2),'g','DisplayName','plant')
     subplot 122;plot(x(:,1),x(:,2),'g','DisplayName','plant')
end


%linearized system
xeq=[1 0]';
A=Amatrix(xeq);
B=[0 1]';
C=[1 0];% measure position (1) or velocity (2)
D=0;
plant=ss(A,B,C,D);

%draws trajectories close to the eigenvectors
t=0:.1:6;
x0=[.244 1.5;1.756 -1.5]';
for i=1:2,
    [~,t,x]=lsim(plant,zeros(size(t)),t,x0(:,i)-xeq);
    x=x+ones(length(x),1)*xeq';%x must be corrected to 
    %account for equilibrium position
    subplot 221;plot(t,x(:,1),'g--','DisplayName','lin plant')
    subplot 223;plot(t,x(:,2),'g--','DisplayName','lin plant')
    subplot 122;plot(x(:,1),x(:,2),'g--','DisplayName','lin plant')
end

%full state, nonlinear
Kfull=[-3.18 0];%[-1.2 -.5];% controller gain
BK=B*Kfull;

%tspan=[0 30];
x0=[.6 1.4;1.5 -1.4]';
for i=1:2,
    [t,x] = ode45(@(t,x) odefun(t,x,BK,xeq), tspan, x0(:,i));
    subplot 221;plot(t,x(:,1),'r','DisplayName','full')
    subplot 223;plot(t,x(:,2),'r','DisplayName','full')
    subplot 122;plot(x(:,1),x(:,2),'r','DisplayName','full')
end

%full state, linear
Afull=A+B*Kfull;
fullstate=ss(Afull,B,C,D);

t=0:.1:30;
%x0=[.96 .01;.97 -.01]';
for i=1:2,
    [~,t,x]=lsim(fullstate,zeros(size(t)),t,x0(:,i)-xeq);
    x=x+ones(length(x),1)*xeq';%x must be corrected to 
    %account for equilibrium position
    subplot 221;plot(t,x(:,1),'r--','DisplayName','lin full')
    subplot 223;plot(t,x(:,2),'r--','DisplayName','lin full')
    subplot 122;plot(x(:,1),x(:,2),'r--','DisplayName','lin full')
end

%observer, nonlinear
BK=B*Kfull;
L=[-3 -4]';%[-.5 -6.17]';%[-.1 -2]';% observer gain

%tspan=[0 30];
x0t=[x0(:,1)' 1*x0(:,1)'-[1 0];x0(:,2)' 1*x0(:,2)'-[1 0]]';
x0t
for i=1:2,
    [t,x] = ode45(@(t,x) odefun2(t,x,BK,L*C,xeq), tspan, x0t(:,i));
    subplot 221;plot(t,x(:,1)-x(:,3),'c','DisplayName','estimate','linewidth',2)
    subplot 223;plot(t,x(:,2),'c','DisplayName','estimate','linewidth',2)
    subplot 122;plot(x(:,1)-x(:,3),x(:,2)-x(:,4),'c','DisplayName','estimate','linewidth',2)

    subplot 221;plot(t,x(:,1),'b','DisplayName','observer')
    subplot 223;plot(t,x(:,2),'b','DisplayName','observer')
    subplot 122;plot(x(:,1),x(:,2),'b','DisplayName','observer')
end

%observer, linear
Ahat=[A+BK -BK;zeros(2) A+L*C];
Bhat=[0 1 0 1]';
Chat=[1 0 0 0];%position, velocity, est. pos., est. vel.
observer=ss(Ahat,Bhat,Chat,D);

t=0:.1:30;
%x0t=[-1 1.327 0.1 0.1;3 -1.326 0.1 0.1]';
for i=1:2,%set to 0 to not plot
    [~,t,x]=lsim(observer,zeros(size(t)),t,x0t(:,i)-[xeq' 0 0]');
    x=x+ones(length(x),1)*[xeq' 0 0];%x must be corrected to 
    %account for equilibrium position
    subplot 221;plot(t,x(:,1)-x(:,3),'c--','DisplayName','lin estimate','linewidth',2)
    subplot 223;plot(t,x(:,2),'c--','DisplayName','lin estimate','linewidth',2)
    subplot 122;plot(x(:,1)-x(:,3),x(:,2)-x(:,4),'c--','DisplayName','lin estimate','linewidth',2)

    subplot 221;plot(t,x(:,1),'b--','DisplayName','lin observer')
    subplot 223;plot(t,x(:,2),'b--','DisplayName','lin observer')
    subplot 122;plot(x(:,1),x(:,2),'b--','DisplayName','lin observer')
end

%final touches on figure
subplot 221;hold off;legend('show');grid;ylabel('position (theta/pi)')
subplot 223;hold off;legend('show');grid;ylabel('velocity (rad/s)');xlabel('time (s)')
subplot 122;hold off;legend('show');grid;ylabel('velocity (rad/s)');xlabel('position(theta/pi)')
xlim([-.1 2.1])
ylim([-1.5 1.5])

%pzmap
figure(2);clf;pzmap(plant,'g',observer,'b',fullstate,'r');grid;legend('show')
end

function fx = f(x)
w02=1;alpha=.2;
fx = [x(2); -w02*sin(pi*x(1))-2*alpha*x(2)];
end

function A = Amatrix(xeq)
w02=1;alpha=.2;
A=[0 1;-w02*cos(pi*xeq(1))*pi -2*alpha];
end

function dxdt = odefun(~,x,BK,xeq)
dxdt = f(x)+BK*(x-xeq);
end

function dxdt = odefun2(~,xt,BK,LC,xeq)
x=xt(1:2);xe=xt(3:4);
dxdt = [f(x);f(x)-f(x-xe)]+[BK -BK;zeros(2) LC]*[x-xeq;xe];
%+[eye(2) zero(2);eye(2) L]*[d n];
end