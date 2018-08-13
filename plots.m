%% 1. Growth rate for one crystal since t0, from D1.
figure('Name', 'Growth rate of spherulites');
subplot(1,2,1);
plot(t_15,G_15(:,1), 'r', 'LineWidth',5);
hold on
plot(t_25,G_25(:,1), 'c', 'LineWidth',3.5);
hold on
plot(t_3,G_3(:,1), 'm', 'LineWidth',4);
hold on
plot(t_4,G_4(:,1), 'b', 'LineWidth',3);
hold on
plot(t_6,G_6(:,1), 'g', 'LineWidth',2);
hold on
plot(t_8,G_8(:,1), 'k', 'LineWidth',1.5);
xlabel('Time (second)');
ylabel('Crystal growth rate ({\mu}m s^-1)');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
ylim([0,0.09]);
box on
subplot(1,2,2);
plot(T_D1_15,G_15(:,1), 'r', 'LineWidth',5);
hold on
plot(T_D1_25,G_25(:,1), 'c', 'LineWidth',3.5);
hold on
plot(T_D1_3,G_3(:,1), 'm', 'LineWidth',4);
hold on
plot(T_D1_4,G_4(:,1), 'b', 'LineWidth',3);
hold on
plot(T_D1_6,G_6(:,1), 'g', 'LineWidth',2);
hold on
plot(T_D1_8,G_8(:,1), 'k', 'LineWidth',1.5);
xlabel('Temperature (^{\circ}C)');
ylabel('Crystal growth rate ({\mu}m s^-1)');
set(gca,'xtick',60:10:170);
set(gca,'xdir','reverse');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
ylim([0,0.09]);
box on


%% 2. Radius growth rate for one crystal since t0, from D1
figure('Name', 'Radius growth rate of spherulites');
subplot(1,2,1);
plot(t_15,delta_Rp_15(:,1), 'r', 'LineWidth',4);
hold on
plot(t_25,delta_Rp_25(:,1), 'c', 'LineWidth',4);
hold on
plot(t_3,delta_Rp_3(:,1), 'm', 'LineWidth',4);
hold on
plot(t_4,delta_Rp_4(:,1), 'b', 'LineWidth',4);
hold on
plot(t_6,delta_Rp_6(:,1), 'g', 'LineWidth',4);
hold on
plot(t_8,delta_Rp_8(:,1), 'k', 'LineWidth',4);
xlabel('Time (second)');
ylabel('Radius growth rate ({\mu}m)');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on
subplot(1,2,2);
plot(T_D1_15,delta_Rp_15(:,1), 'r', 'LineWidth',4);
hold on 
plot(T_D1_25,delta_Rp_25(:,1), 'c', 'LineWidth',4);
hold on
plot(T_D1_3,delta_Rp_3(:,1), 'm', 'LineWidth',4);
hold on
plot(T_D1_4,delta_Rp_4(:,1), 'b', 'LineWidth',4);
hold on
plot(T_D1_6,delta_Rp_6(:,1), 'g', 'LineWidth',4);
hold on
plot(T_D1_8,delta_Rp_8(:,1), 'k', 'LineWidth',4);
xlabel('Temperature (^{\circ}C)');
ylabel('Radius growth rate ({\mu}m)');
set(gca,'xtick',60:10:170);
set(gca,'xdir','reverse');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on

clear G_above_4 G_below_4 G_above_15 G_below_15 G_above_25 G_below_25 G_above_3 G_below_3 G_above_6 G_below_6 G_above_8 G_below_8; 
clear G_4 G_15 G_25 G_3 G_6 G_8;

%% 3. Radius trajectory of single crystal nucleated at t0, for section D1.
figure('Name', 'Radius trajectory of single crystal');
subplot(1,2,1);
plot(t_15,Rp_15(:,1,1), 'r', 'LineWidth',5);
hold on
plot(t_25,Rp_25(:,1,1), 'c', 'LineWidth',3.5);
hold on
plot(t_3,Rp_3(:,1,1), 'm', 'LineWidth',4);
hold on
plot(t_4,Rp_4(:,1,1), 'b', 'LineWidth',3);
hold on
plot(t_6,Rp_6(:,1,1), 'g', 'LineWidth',2);
hold on
plot(t_8,Rp_8(:,1,1), 'k', 'LineWidth',1.5);
xlabel('Time (second)');
ylabel('Accumulative radius ({\mu}m)');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on
subplot(1,2,2);
plot(T_D1_15,Rp_15(:,1,1), 'r', 'LineWidth',5);
hold on
plot(T_D1_25,Rp_25(:,1,1), 'c', 'LineWidth',3.5);
hold on
plot(T_D1_3,Rp_3(:,1,1), 'm', 'LineWidth',4);
hold on
plot(T_D1_4,Rp_4(:,1,1), 'b', 'LineWidth',3);
hold on
plot(T_D1_6,Rp_6(:,1,1), 'g', 'LineWidth',2);
hold on
plot(T_D1_8,Rp_8(:,1,1), 'k', 'LineWidth',1.5);
xlabel('Temperature (^{\circ}C)');
ylabel('Accumulative radius ({\mu}m)');
set(gca,'xtick',60:10:170);
set(gca,'xdir','reverse');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on

clear delta_Rp_4 delta_Rp_15 delta_Rp_25 delta_Rp_3 delta_Rp_6 delta_Rp_8;

%% 4. Volume trajectory of single crystal nucleated at t0, for section D1.
figure('Name', 'Volume trajectory of single crystal');
subplot(1,2,1);
plot(t_15,Vp_15(:,1,1), 'r', 'LineWidth',5);
hold on
plot(t_25,Vp_25(:,1,1), 'c', 'LineWidth',3.5);
hold on
plot(t_3,Vp_3(:,1,1), 'm', 'LineWidth',4);
hold on
plot(t_4,Vp_4(:,1,1), 'b', 'LineWidth',3);
hold on
plot(t_6,Vp_6(:,1,1), 'g', 'LineWidth',2);
hold on
plot(t_8,Vp_8(:,1,1), 'k', 'LineWidth',1.5);
xlabel('Time (second)');
ylabel('Accumulative volume ({\mu}m^3)');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on
subplot(1,2,2);
plot(T_D1_15,Vp_15(:,1,1), 'r', 'LineWidth',5);
hold on
plot(T_D1_25,Vp_25(:,1,1), 'c', 'LineWidth',3.5);
hold on
plot(T_D1_3,Vp_3(:,1,1), 'm', 'LineWidth',4);
hold on
plot(T_D1_4,Vp_4(:,1,1), 'b', 'LineWidth',3);
hold on
plot(T_D1_6,Vp_6(:,1,1), 'g', 'LineWidth',2);
hold on
plot(T_D1_8,Vp_8(:,1,1), 'k', 'LineWidth',1.5);
xlabel('Temperature (^{\circ}C)');
ylabel('Accumulative volume ({\mu}m^3)');
set(gca,'xtick',60:10:170);
set(gca,'xdir','reverse');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on

clear Rp_4 Rp_15 Rp_25 Rp_3 Rp_6 Rp_8;

%% 5. Nucleation numbers in section D1 in each time step
figure('Name', 'The number of nuclei');
subplot(1,2,1);
plot(t_15,N_15(:,1), 'r', t_25,N_25(:,1), 'c', t_3,N_3(:,1), 'm', t_4,N_4(:,1), 'b',t_6,N_6(:,1), 'g', t_8,N_8(:,1), 'k', 'LineWidth',3);
xlabel('Time (second)');
ylabel('The number of nuclei (N)');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on
subplot(1,2,2);
plot(T_15(:,1),N_15(:,1), 'r', T_25(:,1),N_25(:,1), 'c', T_3(:,1),N_3(:,1), 'm', T_4(:,1),N_4(:,1), 'b', T_6(:,1),N_6(:,1), 'g',  T_8(:,1),N_8(:,1), 'k', 'LineWidth',3);
xlabel('Temperature (^{\circ}C)');
ylabel('The number of nuclei (N)');
set(gca,'xtick',60:10:170);set(gca,'xdir','reverse');
legend('0.15 mm', '0.25 mm', '0.3 mm', '0.4 mm', '0.6 mm', '0.8 mm');
box on

%% 6. Total crystal volume in each section
figure('name', 'Crystal volume in each section')
hold on
for tn = 1:667
    y = 0.15 * tn;
    scatter(y, (Vtot_nt_15(end,tn)), 'r', 'LineWidth', 1.5);
end
hold on
for tn = 1:400
    y = 0.25 * tn;
    scatter(y, (Vtot_nt_25(end,tn)), 'c', 'LineWidth', 1.5);
end
hold on
for tn = 1:333
    y = 0.3 * tn;
    scatter(y, (Vtot_nt_3(end,tn)), 'm', 'LineWidth', 1.5);
end
hold on
for tn = 1:250
    y = 0.4 * tn;
    scatter(y, (Vtot_nt_4(end,tn)), 'b', 'LineWidth', 1.5);
end
hold on
for tn = 1:167
    y = 0.6 * tn;
    scatter(y, (Vtot_nt_6(end,tn)), 'g', 'LineWidth',1.5);
end
hold on
for tn = 1:125
    y = 0.8 * tn;
    scatter(y, (Vtot_nt_8(end,tn)), 'k', 'LineWidth', 1.5);
end
hold on
xlabel('Distance from the initial printing point (mm)');
ylabel('Crystal volume ({\mu}m^3)');
box on

%% 7. Total crystal volume trajectory from D1 to the end
figure('name', 'Crystal volume over all sections')
hold on
for tn = 1:667
    V_sum_15 = cumsum(Vtot_nt_15(end, 1:tn));
    y = 0.15 * tn;
    scatter(y, (V_sum_15(end)), 'r', 'LineWidth', 1.5);
end
hold on
for tn = 1:400
    V_sum_25 = cumsum(Vtot_nt_25(end, 1:tn));
    y = 0.25 * tn;
    scatter(y, (V_sum_25(end)), 'c', 'LineWidth', 1.5);
end
hold on
for tn = 1:333
    V_sum_3 = cumsum(Vtot_nt_3(end, 1:tn));
    y = 0.3 * tn;
    scatter(y, (V_sum_3(end)), 'm', 'LineWidth', 1.5);
end
hold on
for tn = 1:250
    V_sum_4 = cumsum(Vtot_nt_4(end, 1:tn));
    y = 0.4 * tn;
    scatter(y, (V_sum_4(end)), 'b', 'LineWidth', 1.5);
end
hold on
for tn = 1:167
    V_sum_6 = cumsum(Vtot_nt_6(end, 1:tn));
    y = 0.6 * tn;
    scatter(y, (V_sum_6(end)), 'g', 'LineWidth', 1.5);
end
hold on
for tn = 1:125
    V_sum_8 = cumsum(Vtot_nt_8(end, 1:tn));
    y = 0.8 * tn;
    scatter(y, (V_sum_8(end)), 'k', 'LineWidth', 1.5);
end
xlabel('Distance from the initial printing point (mm)');
ylabel('Accumulated crystal volume ({\mu}m^3)');
box on