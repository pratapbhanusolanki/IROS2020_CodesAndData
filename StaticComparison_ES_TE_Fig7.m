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
x_hs(1,:) = [17,14];
y_hs(1) = g(acosd((cosd(x_hs(1,1)-xopt1)*cosd(x_hs(1,2)-xopt2))));
k =1;
tf = 10;
t= 0:dt:tf;
xm1 = x_hs(1,1);
xm2 = x_hs(1,2);
x1(k) = xm1 + A*sin(omega*t(k) + phase1);
x2(k) = xm2 + A*sin(omega*t(k) + phase2);

y(k) = g(acosd((cosd(x1(k)-xopt1)*cosd(x2(k)-xopt2))));

theta(1) = randi([-179,180]);
x_hs(2,:) = x_hs(1,:) + r*[cosd(theta(1)),sind(theta(1))];
y_hs(2) = g(acosd((cosd(x_hs(2,1)-xopt1)*cosd(x_hs(2,2)-xopt2))));

del_theta(2) = 120;
theta(2) = theta(1) + del_theta(2);
x_hs(3,:) = x_hs(2,:) + r*[cosd(theta(2)),sind(theta(2))];      %First equilateral triangle
y_hs(3) = g(acosd((cosd(x_hs(3,1)-xopt1)*cosd(x_hs(3,2)-xopt2))));

k = 3;
del_y = (y_hs(k)+y_hs(k-1))/2 - y_hs(k-2);
sgn_y = sign(del_y) + xor(sign(del_y),1); %XOR is used to keep non-zero sign
del_theta(k) = -sgn_y*del_theta(k-1);
theta(k) = theta(k-1) + del_theta(k);



butterorder = 1;
[fb,fa] = butter(butterorder, butterfreq, 'high' );


effort_hs = sum(sum(abs(x_hs)));
effort = 0;

my_line_width = 2;

%For function visualization
C = 5;
step_size = 0.2;
%limit = 20;
x1_contour = -limit:step_size:limit;
x2_contour = -limit:step_size:limit;

for i=1:length(x1_contour)
    for j=1:length(x2_contour)
        Y1(j,i) = g(acosd((cosd(x1_contour(i)-xopt1)*cosd(x2_contour(j)-xopt2))));
    end
end

for k=2:length(t)
    x1(k) = xm1 + A*sin(omega*t(k) + phase1);
    x2(k) = xm2 + A*sin(omega*t(k) + phase2);
    if k>1
        effort = effort + abs(x1(k)-x1(k-1)) + abs(x2(k)-x2(k-1));
    end
    y(k) = g(acosd((cosd(x1(k)-xopt1)*cosd(x2(k)-xopt2))));
    
    lim = max(1,k-10);
    
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
        y_hs(k) = g(acosd((cosd(x_hs(k,1)-xopt1)*cosd(x_hs(k,2)-xopt2))));
        del_y = (y_hs(k)+y_hs(k-1))/2 - y_hs(k-2);
        sgn_y = sign(del_y) + xor(sign(del_y),1); %XOR is used to keep non-zero sign
        del_theta(k) = -sgn_y*del_theta(k-1);
        theta(k) = theta(k-1) + del_theta(k);
        current_effort = abs(r*[cosd(theta(k-1)),sind(theta(k-1))]);
        effort_hs = effort_hs + sum(current_effort);
    end
    
   
end

figure;
h1 = plot(x1,x2,'b--','LineWidth',my_line_width);
hold on;
h2 = plot(x_hs(:,1),x_hs(:,2),'r-','LineWidth',my_line_width);
[dummy,y1_contour] = contour(x1_contour,x2_contour,Y1,'linewidth', 1, 'linecolor','k');
legend([h1 h2], 'Extremum Seeking','Triangular Exploration','Location','NorthWest');
    %drawnow;
axis equal;
xlabel("$$x_1 (^\circ)$$",'interpreter','latex');
ylabel("$$x_2 (^\circ)$$",'interpreter','latex');
