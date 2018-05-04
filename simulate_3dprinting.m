% Parameters
U    = 1500;
Rg   = 1.987;
Goi  = 3.98e7;
Goii = 4.81e11;
Kgi  = 2.55e5;
Kgii = 5.51e5;
No   = 63;
k    = 1 / 54978;
KtoC = 273.15;
Tm   = 212 + KtoC;

% Variables
t_step     = 0.02;         % Time step
t          = 0:t_step:5;   % Time
n_tsteps   = numel(t);     % Number of time steps
n_sections = n_tsteps;     % Number of sections of material deposited
T_D1       = 170 - 22 * t; % Temperature time series for section D1


%% Pre-allocate temperature and mask matrices
% Columns for D1 -- D251, and rows for time.
T = zeros(n_tsteps);                    % Temperature of current section (D1 ... D251) at current time
mask = zeros(n_tsteps, 'logical');      % Mask array to filter data / non-data

for i = 1:n_sections
	last_section_step = n_sections -i +1;
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
% delta_Rp: Change in radius of one crystal at each time step
delta_Rp = G * t_step;

% Rp: Current radius of one particle (crystal) nucleated at t0
Rp_sincet0 = cumsum(delta_Rp, 1, 'omitnan');

% delta_Vp: Change in volume of one particle (crystal) at each time step
delta_Vp = 4/3 * pi * delta_Rp.^3;

%% PLA nucleation
% Total number of particles: not number of nucleated particles at current
% time step
N =  No * exp(-(1/k) ./ ((T + KtoC) .* (Tm - T - KtoC)));

% Number should be zero where T = 0. 
% This is becase where T = 0, it should actually equal NaN instead,
% but this yields an error in the calculation of G.
% Hence use mask afterwards to remove incorrect values.
N(~mask) = 0; 

%% Crystal volume
delta_V = N .* delta_Vp; % Change in total volume of all crystals
Vp_sincet0 = cumsum(delta_Vp, 1, 'omitnan');

%% Total volume accumulated since t = tn.
V = zeros(n_tsteps);
for col = 1:n_sections
	for row = col:n_tsteps
		V(row, col) = sum(delta_Vp(row:end, col), 1);
	end
end

%% Total accumulated volume since t0
Vsec = sum(V, 1); % Total accumulated volume in each section
Vtot = sum(Vsec); % Total accumulated volume in all sections


%% Growth rate for one particle since t0, from D1.
figure('Name', 'Growth rate for D1');
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
plot(t,Rp_sincet0(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius of one particle ({\mu}m)');

%% Cumulative size of one particle since t0, for section D1.
figure('Name', 'Cumulative particle size since t0 for section D1');
subplot(1,2,1);
plot(t,Rp_sincet0(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Radius of one particle ({\mu}m)');
set(gca,'ytick',0:0.01:0.18);
subplot(1,2,2);
plot(t,Vp_sincet0(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Spherical volume of one particle ({\mu}m^3)');

%% Cumulative size of one particle since tn, for section D1.
figure('Name', 'Cumulative particle size since tn for D1');
hold on
for i = 1:n_tsteps
	plot(t(i:end), cumsum(delta_Rp(i:end,1)),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Radius ({\mu}m)');

%% Cumulative volume of one particle since tn, for section D1.
figure('Name', 'Cumulative particle volume since tn for D1');
hold on
for i=1:n_tsteps
	plot(t(i:end), cumsum(delta_Vp(i:end,1)),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
%% Number of particles
figure('Name', 'Number of particles');
subplot(1,2,1);plot(t,N(:,1),'LineWidth',3);
xlabel('Time (second)');
ylabel('Number of spherulities');
axis([0,5,0,25]);set(gca,'xtick',0:0.5:5);set(gca,'ytick',0:2:24);
subplot(1,2,2);plot(T(:,1),N(:,1),'LineWidth',3);
xlabel('Temperature (20^{\circ}C)');
ylabel('Number of spherulities');
axis([60,170,0,25]);set(gca,'xtick',60:10:170);set(gca,'xdir','reverse');set(gca,'ytick',0:2:24)

%% Cumulative volume of all particles since tn, for section D1.
figure('Name', 'Cumulative volume since tn');
hold on
for i=1:n_tsteps
	plot(t(i:end),cumsum(delta_V(i:end,1)),'LineWidth',1.25);
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');

%% Selected cumulative volume since tn for section D1.
figure('Name', 'Selected cumulative volume since tn');
start_times = [1, 26, 51, 101, 151, 201];
n_traces = length(start_times);
hold on
for i = 1:n_traces
	start = start_times(i);
	plot(t(start:end),cumsum(delta_V(start:end,1)), 'LineWidth',3)
end
xlabel('Time (second)');
ylabel('Spherical volume ({\mu}m^3)');
legend('from 0s --> 5s','from 0.5s --> 5s','from 1s --> 5s','from 2s --> 5s','from 3s --> 5s','from 4s --> 5s')	
