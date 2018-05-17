% Parameters
U    = 1500;
Rg   = 1.987;
Goi  = 3.98e7;
Goii = 4.81e11;
Kgi  = 2.55e5;
Kgii = 5.51e5;
No   = 63;
k    = 1 / 54978;
CtoK = 273.15;
Tm   = 212 + CtoK;

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

G_above = Goi * exp(-U ./ (Rg * (T - 30))) .* exp(-Kgi .* (T + CtoK + Tm) ./ (2 * (T + CtoK).^2 .* (Tm - T - CtoK)));
G_below = Goii * exp(-U ./ (Rg * (T - 30))) .* exp(-Kgii .* (T + CtoK + Tm) ./ (2 * (T + CtoK).^2 .* (Tm - T - CtoK)));

G(above_threshold) = G_above(above_threshold);
G(~above_threshold) = G_below(~above_threshold);
G = G / 60;
G(~mask) = NaN;

%% Crystal size
% delta_Rp: Change (increase) in radius of one crystal at each time step
% This is the same for all sizes of crystals (i.e. existing crystals, and
% those nucleated at the current time step).
delta_Rp = G * t_step;

% Rp: Current radius of one particle (crystal)
% Nucleation time (tn) is the page number i.e. Rp(t, D, tn) is the radius
% of one particle in section D, at time t, nucleated at time tn.
Rp = zeros(n_tsteps, n_sections, n_tsteps);
for tn = 1:n_tsteps
    Rp(tn:end, :, tn) = cumsum(delta_Rp(tn:end, :), 1, 'omitnan');
end

% Vp: Volume of one particle (crystal) at each time step
Vp = 4/3 * pi * Rp.^3;

%% PLA nucleation
% Total number of particles: not number of nucleated particles at current
% time step
N =  No * exp(-(1/k) ./ ((T + CtoK) .* (Tm - T - CtoK)));

% Number should be zero where T = 0. 
% This is becase where T = 0, it should actually equal NaN instead,
% but this yields an error in the calculation of G.
% Hence use mask afterwards to remove incorrect values.
N(~mask) = 0; 

%% Accumulated volume
V = zeros(n_tsteps, n_sections, n_tsteps);
% V: Total volume of all crystals
% V(t, D, tn) is the volume of all crystals nucleated at time tn, in
% section D1, at time t.
for tn = 1:n_tsteps
    V(:, :, tn) = N(tn, :) .* Vp(:, :, tn);
end

%% Total volume accumulated since t = 0.
% Vtot: total volume of all particles
Vtot_nt = sum(V, 3);        % Sum of volume over all nucleation times
Vtot_sec = sum(Vtot_nt, 2); % Sum of volume over all sections and nucleation times
Vtot_t = sum(Vtot_nt, 1);   % Sum of volume over all time steps and nucleation times
Vtot_all = sum(Vtot_t);     % Volume since t=0 for all particles.
                            % This is the same as sum(Vtot_sec)

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
