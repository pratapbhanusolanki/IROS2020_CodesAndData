clear all;
close all;
clc;

%For function visualization
global a;
a = 20;
b =  20;

%For function visualization

limit = 25;
num_iteration = 1000 
fs = 100;
f_perturb = fs/3;
omega = f_perturb*pi*2
dt = 1/fs;
phase1 = 0;
phase2 = pi/2;
K= 500;
A = 1;
step_size = sqrt(3)*A;
r = step_size;
butterfreq = 2*(f_perturb/10)/fs;

Mopt = 5;

t= 0:dt:20;
frequencies = [0.001, 0.002, 0.005, 0.01,0.02,0.05, 0.1,0.2,0.5 1,2, 5, 10];
x0 = lhsdesign(1000,2)*2*limit - limit;
for m = 1:length(frequencies)
    f_optima = frequencies(m);
    omega_optima = 2*pi*f_optima;
    m
    for l=398%1:1000
        fprintf('l = %03d', l);
        k =1;
        x_optimal1(k) = Mopt*r*cos(omega*t(k));
        x_optimal2(k) = Mopt*r*sin(omega*t(k));
        %x_hs(1,:) = [randi([-limit,limit]),randi([-limit,limit])];
        x_hs(1,:) = x0(l,:);
        y_hs(1) = g(acosd((cosd(x_hs(1,1)-x_optimal1(k))*cosd(x_hs(1,2)-x_optimal2(k)))));

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

        
        for k=1:length(t)
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
            dist_tri(k) = norm(x_hs(k,:)-[x_optimal1(k),x_optimal2(k)]);
            dist_ext(k) = norm([x1(k)-x_optimal1(k),x2(k)-x_optimal2(k)]);
        end
        fprintf('\b\b\b\b\b\b\b');
        y_ext(l) = mean(y(k+1-100:k));
        dist_mean_ext(l) = mean(dist_ext(k+1-100:k));
        
        y_tri(l) = mean(y_hs(k+1-100:k));
        dist_mean_tri(l) = mean(dist_tri(k+1-100:k));
        
    end
    y_ext_f_mean(m) = mean(y_ext);
    y_ext_f_std(m) = std(y_ext);
    
    dist_ext_f_mean(m) = mean(dist_mean_ext);
    dist_ext_f_std(m) = std(dist_mean_ext);
    
    y_tri_f_mean(m) = mean(y_tri);
    y_tri_f_std(m) = std(y_tri);
    
    dist_tri_f_mean(m) = mean(dist_mean_tri);
    dist_tri_f_std(m) = std(dist_mean_tri);
    
end

figure('units','normalized','outerposition',[0 0 1 1]);
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
%subplot(2,1,2)
%export_fig('Media/moving_optima.eps');






figure('units','normalized','outerposition',[0 0 1 1]);
plot(x1,x2);
hold on;
plot(x_hs(:,1),x_hs(:,2),'r--');
plot(x_optimal1,x_optimal2,'k-.');



%For function visualization
C = 5;
step_size = 0.2;
limit = 20;
x1 = -limit:step_size:limit;
x2 = -limit:step_size:limit;

for i=1:length(x1)
    for j=1:length(x2)
        Y1(j,i) = g(acosd((cosd(x1(i))*cosd(x2(j)))));
    end
end

font_size = 24;
gca_font_size = 18;
[dummy,y1_contour] = contour(x1,x2,Y1,'linewidth', 1, 'linecolor','k');
set(gca,'FontSize',gca_font_size);
xlabel('$$x_1$$','Interpreter','latex','FontSize', 24);
ylabel('$$x_2$$','Interpreter','latex','FontSize', 24);
%export_fig('Media/moving_optima.eps','-transparent');
%legend('Triangular Exploration', 'Extremum Seeking')
effort
effort_hs
figure;
h1 = errorbar(frequencies, y_ext_f_mean, y_ext_f_std);
hold on; 
h2= errorbar(frequencies, y_tri_f_mean,y_tri_f_std,'r');
set(get(h1,'Parent'), 'XScale', 'log');
set(get(h2,'Parent'), 'XScale', 'log');
set(gca,'FontSize',gca_font_size);
xlabel('$$f_M$$ (Hz)','Interpreter','latex','FontSize', 24);
ylabel('$$y_{ave}$$ (V)','Interpreter','latex','FontSize', 24);
figure;
h1 = errorbar(frequencies, dist_ext_f_mean, dist_ext_f_std);
hold on; 
h2 = errorbar(frequencies, dist_tri_f_mean, dist_tri_f_std,'r');
set(get(h1,'Parent'), 'XScale', 'log');
set(get(h2,'Parent'), 'XScale', 'log');
xlabel('$$f_M$$ (Hz)','Interpreter','latex','FontSize', 24);
ylabel('$$r_{ave}$$ (degree)','Interpreter','latex','FontSize', 24);

