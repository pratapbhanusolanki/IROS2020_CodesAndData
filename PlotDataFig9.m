clear all;
close all;
clc;

load DataFrequencyResponseFig9.mat
figure;
h1 = errorbar(frequencies, y_ext_f_mean, y_ext_f_std);
hold on; 
h2= errorbar(frequencies, y_tri_f_mean,y_tri_f_std,'r');
set(get(h1,'Parent'), 'XScale', 'log');
set(get(h2,'Parent'), 'XScale', 'log');
set(gca,'FontSize',gca_font_size);
xlabel('$$f_M$$ (Hz)','Interpreter','latex','FontSize', 24);
ylabel('$$y_{ave}$$ (V)','Interpreter','latex','FontSize', 24);
xlim([0.001,20]);
legend('Extremum Seeking','Triangular Exploration' );
figure;
h1 = errorbar(frequencies, dist_ext_f_mean, dist_ext_f_std,'b--');
hold on; 
h2 = errorbar(frequencies, dist_tri_f_mean, dist_tri_f_std,'r-');
set(get(h1,'Parent'), 'XScale', 'log');
set(get(h2,'Parent'), 'XScale', 'log');
xlabel('$$f_M$$ (Hz)','Interpreter','latex','FontSize', 24);
ylabel('$$r_{ave} (^o)$$','Interpreter','latex','FontSize', 24);
ylim([0,80]);
xlim([0.001,20]);
set(gca,'FontSize',gca_font_size);
legend('Extremum Seeking','Triangular Exploration' );

