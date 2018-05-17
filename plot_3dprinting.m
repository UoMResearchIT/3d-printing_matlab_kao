%% Growth rate for one particle since t0, from D1.
figure('Name', 'Growth rate for particles in section D1');
subplot(1,3,1);
plot(t,G(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Growth rate ({\mu}m s^-1)');
subplot(1,3,2);
plot(T_D1,G(:,1),'LineWidth',3);
xlabel('Temperature (20^{\circ}C)');
ylabel('Growth rate ({\mu}m s^-1)');
axis([60,170,0,0.08]);
set(gca,'xtick',60:50:170);
set(gca,'xdir','reverse')
subplot(1,3,3);
plot(t,Rp(:,1,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius of one particle ({\mu}m)');

%% Size trajectory of one particle nucleated at t0, for section D1.
figure('Name', 'Size trajectory for a particle nucleated at t0 in section D1');
subplot(1,2,1);
plot(t,Rp(:,1,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius of one particle ({\mu}m)');
set(gca,'ytick',0:0.01:0.18);
subplot(1,2,2);
plot(t,Vp(:,1,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Spherical volume of one particle ({\mu}m^3)');

%% Size trajectory of one particle nucleated at tn, in section D1.
figure('Name', 'Size trajectory for a particle nucleated at tn in section D1');
hold on
for tn = 1:n_tsteps
	plot(t(tn:end), Rp(tn:end, 1, tn), 'LineWidth', 1.25);
end
xlabel('Time (second)');
ylabel('Radius ({\mu}m)');

%% Volume trajectory of one particle nucleated at tn, in section D1.
figure('Name', 'Volume trajectory of one particle nucleated at tn in section D1');
hold on
for tn=1:n_tsteps
	plot(t(tn:end), Vp(tn:end, 1, tn), 'LineWidth', 1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
%% Number of particles nucleated in section D1 in each time step
figure('Name', 'Number of particles nucleated in section D1');
subplot(1,2,1);plot(t,N(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Number of particles nucleated');
axis([0,5,0,25]);set(gca,'xtick',0:0.5:5);set(gca,'ytick',0:2:24);
subplot(1,2,2);plot(T(:,1),N(:,1),'LineWidth',3);
xlabel('Temperature (20^{\circ}C)');
ylabel('Number of particles nucleated');
axis([60,170,0,25]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse');set(gca,'ytick',0:2:24)

%% Total volume of all particles nucleated at tn, in section D1.
figure('Name', 'Volume of particles nucleated at time tn in section D1');
hold on
for nt=1:n_tsteps
	plot(t(nt:end), V(nt:end, 1, nt), 'LineWidth', 1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');

%% Volume trajectory of all particles nucleated at selected tn in section D1.
figure('Name', 'Volume of particles nucleated at selected times in section D1');
start_times = [1, 26, 51, 101, 151, 201];
n_traces = length(start_times);
hold on
for tn = [1, 26, 51, 101, 151, 201]
	plot(t(tn:end), V(tn:end, 1, tn), 'LineWidth', 3)
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
legend('from 0s --> 5s','from 0.5s --> 5s','from 1s --> 5s','from 2s --> 5s','from 3s --> 5s','from 4s --> 5s')	

%% Animate volume trajectory from each nucleation step
figure('Name', 'Total volume animation')
for tn = 1:n_tsteps
    imagesc(V(:, :, tn))
    title(['Nucleation time = ', num2str(t(tn))])
    xlabel('Section')
    ylabel('Time step')
    c = colorbar;
    c.Label.String = 'Volume of particles nucleated at current time step';
    caxis([0, 0.08])
    drawnow
end