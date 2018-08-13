%% Variables for simulation
t_step_4     = 0.02;         % Time step
t_4          = 0:t_step_4:33;   % Time
n_tsteps_4   = numel(t_4);     % Number of time steps
n_sections_4 = n_tsteps_4;     % Number of sections of material deposited
T_D1_4       = 170 - 3.33 * t_4; % Temperature time series for section D1

%% Pre-allocate temperature and mask matrices
% Columns for D1 -- to end, and rows for time.
T_4 = zeros(n_tsteps_4);                    % Temperature of current section (D1 ... end) at current time
mask_4 = zeros(n_tsteps_4, 'logical');      % Mask array to filter data / non-data
for i = 1:250
	last_section_step_4 = n_sections_4 -i +1;
	T_4(i:end, i) = T_D1_4(1:last_section_step_4);
	mask_4(i:end, i) = true;
end

%% Growth rate
T_threshold_4 = 120;   % Temperature threshold
G_4 = zeros(n_tsteps_4); % Growth rate
above_threshold_4 = T_4 >= T_threshold_4;

G_above_4 = Goi * exp(-U ./ (Rg * (T_4 - 30))) .* exp(-Kgi .* (T_4 + C2K + Tm) ./ (2 * (T_4 + C2K).^2 .* (Tm - T_4 - C2K)));
G_below_4 = Goii * exp(-U ./ (Rg * (T_4 - 30))) .* exp(-Kgii .* (T_4 + C2K + Tm) ./ (2 * (T_4 + C2K).^2 .* (Tm - T_4 - C2K)));

G_4(above_threshold_4) = G_above_4(above_threshold_4);
G_4(~above_threshold_4) = G_below_4(~above_threshold_4);
G_4 = G_4 / 60;
G_4(~mask_4) = NaN;

%% Crystal size
% delta_Rp: Change (increase) in radius of one crystal at each time step
% This is the same for all sizes of crystals (i.e. existing crystals, and
% those nucleated at the current time step).
delta_Rp_4 = G_4 * t_step_4;

%% This is the section where MATLAB runs out of memory
% Rp: Current radius of one particle (crystal)
% Nucleation time (tn) is the page number i.e. Rp(t, D, tn) is the radius
% of one particle in section D, at time t, nucleated at time tn.
Rp_4 = zeros(n_tsteps_4, n_sections_4, n_tsteps_4);
for tn = 1:n_tsteps_4
    Rp_4(tn:end, :, tn) = cumsum(delta_Rp_4(tn:end, :), 1, 'omitnan');
end

% Vp: Volume of one particle (crystal) at each time step
Vp_4 = 4/3 * pi * Rp_4.^3;

%% PLA nucleation
% Total number of particles: not number of nucleated particles at current
% time step
N_4 =  No_4 * exp(-(1/k) ./ ((T_4 + C2K) .* (Tm - T_4 - C2K)));

% Number should be zero where T = 0. 
% This is becase where T = 0, it should actually equal NaN instead,
% but this yields an error in the calculation of G.
% Hence use mask afterwards to remove incorrect values.
N_4(~mask_4) = 0;

%% Accumulated volume
% V: Total volume of all crystals
% Vtot_nt(t, D) is the volume of all crystals summed over all nucleation
% times (tn) in section D1, at time t.
Vtot_nt_4 = zeros(n_tsteps_4, n_sections_4);
for tn = 1:n_tsteps_4
    Vtot_nt_4 = Vtot_nt_4 + N_4(tn, :) .* Vp_4(:, :, tn)
end

%% Total volume accumulated since t = 0.
% Vtot: total volume of all particles
Vtot_sec_4 = sum(Vtot_nt_4, 2); % Sum of volume over all sections and nucleation times
Vtot_t_4 = sum(Vtot_nt_4, 1);   % Sum of volume over all time steps and nucleation times
Vtot_all_4 = sum(Vtot_t_4);     % Volume since t=0 for all particles.
                                % This is the same as sum(Vtot_sec)