clc; clear all; close all;

%% Balanç calor
%balanc_calor = readmatrix('estudi tipus fisic/balanc_calor_T0_21.csv','NumHeaderLines',1);
balanc_calor_roure = readmatrix('estudi tipus fisic/balanc_calor_roure.csv','NumHeaderLines',1);
balanc_calor_roure_15 = readmatrix('estudi tipus fisic/balanc_calor_roure_15.csv','NumHeaderLines',1);
balanc_calor_roure_2 = readmatrix('estudi tipus fisic/balanc_calor_roure_2.csv','NumHeaderLines',1);
balanc_calor_roure_3 = readmatrix('estudi tipus fisic/balanc_calor_roure_3.csv','NumHeaderLines',1);
balanc_calor_roure_4 = readmatrix('estudi tipus fisic/balanc_calor_roure_4.csv','NumHeaderLines',1);
%%
t = balanc_calor(2,:);
Qin = balanc_calor(3,:);
Qout = balanc_calor(4,:);
Qacc = balanc_calor(5,:);
Qtotal = balanc_calor(6,:);
Qin_roure = balanc_calor_roure(3,:);
Qout_roure = balanc_calor_roure(4,:);
Qacc_roure = balanc_calor_roure(5,:);
Qtotal_roure = balanc_calor_roure(6,:);
Qout_roure_15 = balanc_calor_roure_15(4,:);
Qout_roure_2 = balanc_calor_roure_2(4,:);
Qout_roure_3 = balanc_calor_roure_3(4,:);
Qout_roure_4 = balanc_calor_roure_4(4,:);
%%
max(abs(Qtotal))

%% Balanç global
plot(t,Qtotal,'b')
title('Balan\c{c} global d''energia en cada temps','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Energia per unitat de temps [W]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
%legend('','Interpreter','latex', 'Location','best');
grid minor, grid on;
%% Energia intercanviada x=0
plot(t,Qin,'g')
title('Energia intercanviada en $x = 0$','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Flux d''energia per unitat de temps [W]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
%legend('','Interpreter','latex', 'Location','best');
grid minor, grid on;
%%
purple = [139 0 255]/255;
hold on;
plot(t,Qout_roure_4,'r',t,Qout_roure_3,'g',t,Qout_roure_2,'c',t,Qout_roure_15,'b')
plot(t,Qout_roure,'color',purple)
hold off;
legend('$e=4\,$m','$e=3\,$m','$e=2\,$m','$e=1.5\,$m','$e=1\,$m','Interpreter','latex', 'Location','southeast');
title('Energia intercanviada en $x = e$','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Flux d''energia per unitat de temps [W]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
grid minor, grid on;

%% Energia intercanviada x=0 i x=e zoom
plot(t,Qin_roure,'g')
hold on;
plot(t,Qout_roure,'r')
%xline(3600*24,'--');
%xlim([0 1e5])
title('Energia intercanviada','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Flux d''energia per unitat de temps [W]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('$x=0$','$x=e$','Interpreter','latex', 'Location','best');
grid minor, grid on;
%%
min(Qout_roure)
min(Qout)
%%
%T_numerica = readmatrix('estudi tipus fisic/resultats_temperatura_T0_21.csv','NumHeaderLines',1);
T_numerica_roure = readmatrix('estudi tipus fisic/resultats_temperatura_roure_2.csv','NumHeaderLines',1);
x = T_numerica_roure(1:end,1);
T_numerica = T_numerica(1:end,2:end);
T_numerica_roure = T_numerica_roure(1:end,2:end);
%%
plot(t,T_numerica_roure(1,:),'-g',t,T_numerica_roure(27,:),'-b',t,T_numerica_roure(end,:),'-r')
%xline(3600*24,'--');
%xlim([0 1e5])
title('Evoluci\''o de la temperatura','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Temperatura [$^\circ$C]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend(sprintf('x = %d m',x(1)),sprintf('x = %.0f m',x(27)),sprintf('$x = %d$ m',x(end)),'Interpreter','latex', 'Location','best');
grid minor, grid on;