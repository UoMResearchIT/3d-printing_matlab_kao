% Parameters
U = 1500;
Rg = 1.987;
Goi = 3.98e7;
Goii = 4.81e11;
Kgi = 2.55e5;
Kgii = 5.51e5;
Tm = 212 + 273.15;
No = 63;
k = 1 / 54978;
KtoC = 273.15;

t_step = 0.02;
t = 0:t_step:5;
n_tsteps = numel(t);
T_D1 = 170 - 22 * t; % Temperature time series for section D1


%% Create matrices with columns for D1 -- D251, and rows for time.
% Pre-allocate with zeros.
delta_Rp = zeros(n_tsteps);  % Change in radius of one crystal at each time step
Rp = zeros(n_tsteps);        % Current radius of one particle (crystal) at current time step
Vp = zeros(n_tsteps);        % Volume of one particle (crystal) at current time step
T = zeros(n_tsteps);         % Temperature of current section (D1 ... D251) at current time
mask = zeros(n_tsteps, 'logical');      % Mask array to filter data / non-data

for i = 1:n_tsteps
	last_section_step = n_tsteps -i +1;
	T(i:end, i) = T_D1(1:last_section_step);
	mask(i:end, i) = true;
end

%% Growth rate
T_threshold = 120;   % Temperature threshold
G = zeros(n_tsteps); % Growth rate
above_threshold = T >= T_threshold;

G_above = Goi * exp(-U ./ (Rg * (T - 30))) .* exp(-Kgi .* (T + KtoC + Tm) ./ (2 * (T + KtoC).^2 .* (Tm - T - KtoC)));
G_below = Goii * exp(-U ./ (Rg * (T - 30))) .* exp(-Kgii .* (T + KtoC + Tm) ./ (2 * (T + KtoC).^2 .* (Tm - T - KtoC)));

G(above_threshold) = G_above(above_threshold);
G(~above_threshold) = G_below(~above_threshold);
G = G / 60;
G(~mask) = NaN;

%% Crystal size
delta_Rp = G * t_step;
Rp = cumsum(delta_Rp, 1, 'omitnan');
Vp = 4/3 * pi * Rp.^3;

%% PLA nucleation
N =  No * exp(-(1/k) ./ ((T + KtoC) .* (Tm - T - KtoC)));
% Number should be zero where T = 0. 
% Where T = 0, T should really equal NaN instead, but this yields an error
% in the calculation of G.
N(~mask) = 0; 

%% Crystal volume
V = N .* 4/3 * pi * Rp .^3; % Total volume of all crystals

%% Cumulative sums
% Calculate cumulative sums of radius of one particle, and total crystal
% volume.
% This doesn't make sense to me -- are you sure this is what you want?

Rp_sum = cumsum(Rp, 1);
V_sum = cumsum(V, 1);



%% Plotting
figure(1);
subplot(1,3,1);plot(t,G(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Growth rate ({\mu}m s^-1)');
subplot(1,3,2);plot(T_D1,G(:,1),'LineWidth',3);
xlabel('Temperature (20^{\circ}C)');
ylabel('Growth rate ({\mu}m s^-1)');
axis([60,170,0,0.08]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse')
subplot(1,3,3);plot(t,Rp(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius ({\mu}m)');

figure(2);
subplot(1,2,1);plot(t,Rp(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius ({\mu}m)');
set(gca,'ytick',0:0.01:0.18);
subplot(1,2,2);plot(t,Vp(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
set(gca,'ytick',0:0.001:0.025);

figure(3);
hold on
for i = 1:n_tsteps
	plot(t(i:end), Rp(i:end,i),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Radius ({\mu}m)');
%set(gca,'ytick',0:0.01:0.18);

figure(4);
hold on
for i=1:n_tsteps
	plot(t(i:end), V(i:end,i),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
set(gca,'ytick',0:0.002:0.024);

figure(5);
subplot(1,2,1);plot(t,N(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Number of spherulities');
axis([0,5,0,25]);set(gca,'xtick',0:0.5:5);set(gca,'ytick',0:2:24);
subplot(1,2,2);plot(T(:,1),N(:,1),'LineWidth',3);
xlabel('Temperature (20^{\circ}C)');
ylabel('Number of spherulities');
axis([60,170,0,25]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse');set(gca,'ytick',0:2:24)
%%
figure(6);
hold on
for i=1:n_tsteps
	plot(t(i:end),V_sum(i:end,i),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
set(gca,'ytick',0:0.01:0.15);
%%
figure(7);
start_times = [1, 26, 51, 101, 201];
n_traces = length(start_times);
hold on
for i = 1:n_traces
	start = start_times(i);
	plot(t(start:end),V_sum(start:end,start), 'LineWidth',3)
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
legend('from 0s --> 5s','from 0.5s --> 5s','from 1s --> 5s','from 2s --> 5s','from 3s --> 5s','from 4s --> 5s')	
