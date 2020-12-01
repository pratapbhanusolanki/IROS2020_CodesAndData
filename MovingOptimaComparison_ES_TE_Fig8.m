clear all;
close all;
clc;

%For function visualization
global a;
a = 30;
b =  20;

%For function visualization

limit = 20;
num_iteration = 1000 
fs = 100;
f_perturb = fs/3;
omega = f_perturb*pi*2;
dt = 1/fs;
phase1 = 0;
phase2 = pi/2;
K= 500;
A = 1;
step_size = sqrt(3)*A;
r = step_size;
butterfreq = 2*(f_perturb/10)/fs;

xopt1 = 4;
xopt2 = -7;

%x_hs(1,:) = [randi([-limit,limit]),randi([-limit,limit])];
x_hs(1,:) = [17,-14];
k =1;
tf = 20;
t= 0:dt:tf;
k =1;
Mopt = 5;
x_optimal1(k) = Mopt*r*cos(omega*t(k));
x_optimal2(k) = Mopt*r*sin(omega*t(k));
%x_hs(1,:) = [randi([-limit,limit]),randi([-limit,limit])];
%x_hs(1,:) = x0(l,:);
y_hs(1) = g(acosd((cosd(x_hs(1,1)-x_optimal1(k))*cosd(x_hs(1,2)-x_optimal2(k)))));

xm1 = x_hs(1,1);
xm2 = x_hs(1,2);
x1(k) = xm1 + A*sin(omega*t(k) + phase1);
x2(k) = xm2 + A*sin(omega*t(k) + phase2);
y(k) = g(acosd((cosd(x1(k)-x_optimal1(k))*cosd(x2(k)-x_optimal2(k)))));

theta(2) =  randi([0,360]);
k = 2;
x_optimal1(k) = Mopt*r*cos(omega*t(k));
x_optimal2(k) = Mopt*r*sin(omega*t(k));
x_hs(2,:) = x_hs(1,:) + r*[cosd(theta(2)),sind(theta(2))];
y_hs(2) = g(acosd((cosd(x_hs(2,1)-x_optimal1(k))*cosd(x_hs(2,2)-x_optimal2(k)))));

del_theta(3) = 120;
k = 3;
x_optimal1(k) = Mopt*r*cos(omega*t(k));
x_optimal2(k) = Mopt*r*sin(omega*t(k));
theta(3) = theta(2) + del_theta(3);
x_hs(3,:) = x_hs(2,:) + r*[cosd(theta(3)),sind(theta(3))];      %First equilateral triangle
y_hs(3) = g(acosd((cosd(x_hs(3,1)-x_optimal1(k))*cosd(x_hs(3,2)-x_optimal2(k)))));


xm1 = x_hs(1,1);
xm2 = x_hs(1,2);



butterorder = 1;
[fb,fa] = butter(butterorder, butterfreq, 'high' );


effort_hs = sum(sum(abs(x_hs)));
effort = 0;

my_line_width = 2;
f_optima = 0.05;
omega_optima = 2*pi*f_optima;

for k=2:length(t)
    x_optimal1(k) = Mopt*r*cos(omega_optima*t(k));
    x_optimal2(k) = Mopt*r*sin(omega_optima*t(k));

    x1(k) = xm1 + A*sin(omega*t(k) + phase1);
    x2(k) = xm2 + A*sin(omega*t(k) + phase2);
    if k>1
        effort = effort + abs(x1(k)-x1(k-1)) + abs(x2(k)-x2(k-1));
    end
    %y(k) = g(x1(k))*g(x2(k));
    y(k) = g(acosd((cosd(x1(k)-x_optimal1(k))*cosd(x2(k)-x_optimal2(k)))));
    lim = max(1,k-100);

    %yh = filter(fb,fa,y(lim:end));
    yh = y(k) - mean(y(lim:k));
    y_hpf(k) = yh;

    zeta1 = y_hpf(k)*sin(omega*t(k) + phase1);
    zeta2 = y_hpf(k)*sin(omega*t(k) + phase2);

    effort1 =  K*zeta1*dt;
    effort2 =  K*zeta2*dt;
    xm1 = xm1 + effort1;
    xm2 = xm2 + effort2;
    effort = effort + abs(effort1) + abs(effort2);

    if(k>3)
        x_hs(k,:) = x_hs(k-1,:) + r*[cosd(theta(k-1)),sind(theta(k-1))]; 
        y_hs(k) = g(acosd((cosd(x_hs(k,1)-x_optimal1(k))*cosd(x_hs(k,2)-x_optimal2(k)))));
        del_y = (y_hs(k)+y_hs(k-1))/2 - y_hs(k-2);
        sgn_y = sign(del_y) + xor(sign(del_y),1); %XOR is used to keep non-zero sign
        del_theta(k) = -sgn_y*del_theta(k-1);
        theta(k) = theta(k-1) + del_theta(k);
        current_effort = abs(r*[cosd(theta(k-1)),sind(theta(k-1))]);
        effort_hs = effort_hs + sum(current_effort);
    end
end

figure();
subplot(3,1,1);
plot(t,x_hs(:,1),'r--');
hold on;
plot(t,x1,'b-.');
plot(t,x_optimal1,'k');
legend('Triangular Exploration', 'Extremum Seeking','Optimum point')
gca_font_size = 18;
set(gca,'FontSize',gca_font_size);
set(gca,'xticklabel',{[]});
ylabel("$$x_1 (^o)$$",'interpreter','latex');

subplot(3,1,2);
plot(t,x_hs(:,2),'r--');

hold on;
plot(t,x2,'b-.');
plot(t,x_optimal2,'k-');
gca_font_size = 18;
set(gca,'FontSize',gca_font_size);
set(gca,'xticklabel',{[]});
ylabel("$$x_2 (^o)$$",'interpreter','latex');

subplot(3,1,3);
plot(t,y,'b-.');
hold on;
plot(t,y_hs,'r--');
xlabel("$$t$$ (s)",'interpreter','latex');
ylabel("$$y$$ (V)",'interpreter','latex');
ylim([0,1.1]);
font_size = 24;
gca_font_size = 18;
set(gca,'FontSize',gca_font_size);
legend('Triangular Exploration', 'Extremum Seeking')