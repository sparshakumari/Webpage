function mean_power
clc
close all
clear all

%########## CONTROL GROUP ##################

P_c = readmatrix('Control_power_values_10_epochs.xlsx');
F_c = P_c(:,1);
data_c = P_c;
data_c(:,1) = [];
meanpxx_c = mean(data_c,2);                                                     %calculates mean of the rows

sem_c = std(data_c,[],2)/sqrt(length(data_c));                                         % calculates the standard error mean

%########## TEST GROUP ##################

P_t = readmatrix('DGAL_power_values_10_epochs.xlsx');
F_t = P_t(:,1);
data_t = P_t;
data_t(:,1) = [];
meanpxx_t = mean(data_t,2);                                                     %calculates mean of the rows

sem_t = std(data_t,[],2)/sqrt(length(data_t));                                         % calculates the standard error mean

%########## CONTROL GROUP ##################

figure('WindowState','maximized','Color',[1 1 1])
plot(F_c,meanpxx_c, 'k','LineWidth',2.0)                                               % plot power of the individual theta region
set(gca,'box','off')
hold on
e_c = errorbar(F_c, meanpxx_c, sem_c,'LineStyle','none','Color', 'k','LineWidth',1.0, 'HandleVisibility','off')
e_c.Marker = 'o';
e_c.MarkerSize = 5;
e_c.MarkerEdgeColor = 'black';
e_c.MarkerFaceColor = 'white';

%########## TEST GROUP ##################
hold on
plot(F_t,meanpxx_t, 'r','LineWidth',2.0)                                               % plot power of the individual theta region
set(gca,'box','off')
hold on
e_t = errorbar(F_t, meanpxx_t, sem_t,'LineStyle','none','Color', 'k','LineWidth',1.0, 'HandleVisibility','off')
e_t.Marker = 'o';
e_t.MarkerSize = 5;
e_t.MarkerEdgeColor = 'black';
e_t.MarkerFaceColor = 'white';

xlim([0 30])
ylim([0 0.005])
legend('Control','D-Galactose')
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');

end