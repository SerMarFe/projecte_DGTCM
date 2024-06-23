clc; clear all; close all;
%%
T_analitica = readmatrix('resultats_analitics_temperatura.csv');
T_numerica = readmatrix('resultats_temperatura.csv');
x = T_numerica(1:end,1);
T_numerica = T_numerica(1:end,2:end);
%%
set(0, 'DefaultTextInterpreter', 'latex');

%% Solució analítica
t = 0:1:1e5;
plot(x,T_analitica(:,1:end),'-r',LineWidth=1);
title('Soluci\''o anal\''itica vs num\`erica','Interpreter','latex');
xlabel('Posici\''o en la placa [m]');
ylabel('Temperatura [$^\circ$C]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
legend('Anal\''itic','Interpreter','latex', 'Location','best');
grid minor, grid on;

hold on;
newcolors = ["#0B0  " "#00F" "#50F" "#A0F"];
colors = [
    255, 165, 0;   % #ef476f
    255, 209, 102;  % #ffd166
    6, 214, 160;    % #06d6a0
    17, 138, 178;   % #118ab2
    7, 59, 76       % #073b4c
] / 255;

colororder(colors)
plot(x,T_numerica(:,[1,floor(end/32),floor(end/16),end]),'--')

legendText0 = sprintf('t = %d s', 0);
legendText1 = sprintf('t = %d s', floor(length(T_numerica)/32));
legendText2 = sprintf('t = %d s', floor(length(T_numerica)/16));
legendText3 = sprintf('t = %d s', floor(length(T_numerica)));
% Add legend to the plot
legend('Anal\''itic', legendText0, legendText1, legendText2, legendText3);


%% Balanç calor
balanc_calor = readmatrix('balanc_calor.csv','NumHeaderLines',1);
%%
t = balanc_calor(2,:);
Qin = balanc_calor(3,:);
Qout = balanc_calor(4,:);
Qacc = balanc_calor(5,:);
Qtotal = balanc_calor(6,:);
%% Balanç global
plot(t,Qtotal,'b')
title('Balan\c{c} global d''energia en cada temps','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Flux d''energia per unitat de temps [W]');
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
%% Energia intercanviada x=e
plot(t,Qout,'r')
title('Energia intercanviada en $x = e$','Interpreter','latex');
xlabel('Temps [s]');
ylabel('Flux d''energia per unitat de temps [W]');
box on;
set(gca,'TickLabelInterpreter','latex','XLim',xlim*1,'YLim',ylim*1);
%legend('','Interpreter','latex', 'Location','best');
grid minor, grid on;
