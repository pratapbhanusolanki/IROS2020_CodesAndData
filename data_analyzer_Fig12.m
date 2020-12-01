%Uses a simpler model g(x,y) = g(x)*g(y) in the simulation 

clear all;
close all;
clc;
addpath('npy-matlab-master')
unzip('real_time_data_underwater_bot_3D_2020-02-27_01:28_AM.npz','temp_data');
y = readNPY('temp_data/y_all.npy');
u = readNPY('temp_data/u_all.npy');
theta = readNPY('temp_data/theta.npy');
psi = readNPY('temp_data/psi.npy');

x = [7;0];
for i=2:length(y)
    x(:,i) = x(:,i-1) + u(i,2:3)';
end

% figure;
% plot(x(1,:),x(2,:));
% axis equal;

load BitTransmissionDataIROS_Fig12.mat

t = 0.5*(1:length(NumBitsArray));


font_size = 24;
gca_font_size = 18;


figure('units','normalized','outerposition',[0 0 0.5 1]);
subplot(2,1,1);
%plot(t,4*y,'b');
set(gca,'FontSize',gca_font_size);
ylabel('$$y$$ (V)','Interpreter','latex','FontSize', font_size);
xlim([0,t(end)]);
ylim([0,4.2*max(y)]);
hold on;

set(gca,'xticklabel',{[]});
subplot(2,1,2);
%plot(t(1:end-1),2*NumBitsArray(1:end-1),'b-');
%hold on;
%plot(t(1:end-1),2*ErrBitsArray(1:end-1),'r-.');
set(gca,'FontSize',gca_font_size);
xlabel('$$t$$ (s)','Interpreter','latex','FontSize', font_size);
ylabel('Data rates (bits/s)','Interpreter','latex','FontSize', font_size);
xlim([0,t(end)]);
ylim([0,2*max(NumBitsArray)]);
hold on;


%plot(t(1:end-1),2*NumBitsArray(1:end-1)-2*ErrBitsArray(1:end-1),'k--');

%legend("Received bps", "Error bps", "Effective bps");
my_line_width = 2;
for n = 1:length(t)-1
    subplot(2,1,1);
    h1 = plot(t(1:n),4*y(1:n),'b','LineWidth',my_line_width);
    set(gca,'FontSize',gca_font_size);
    set(gca,'xticklabel',{[]});
    
    subplot(2,1,2);
    h1 = plot(t(1:n),2*NumBitsArray(1:n),'b-','LineWidth',my_line_width);
    hold on;
    h2 = plot(t(1:n),2*ErrBitsArray(1:n),'r-.','LineWidth',my_line_width);
    h3 = plot(t(1:n),2*NumBitsArray(1:n)-2*ErrBitsArray(1:n),'k--');
    %legend("Received bps", "Error bps", "Effective bps");
    legend([h1 h2 h3], "Received data-rate", "Error data-rate", "Effective data-rate",'Location','southeast')
    %drawnow;
   % M(n) = getframe(gcf);
end

 
% v = VideoWriter('./experimental_run.mp4');
% v.FrameRate = 2;
% open(v);
% writeVideo(v,M)
% close(v);
