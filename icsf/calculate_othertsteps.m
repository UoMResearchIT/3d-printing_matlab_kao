t_step_15     = 0.0075;         
t_15          = 0:t_step_15:33;   
n_tsteps_15   = numel(t_15);     
n_sections_15 = n_tsteps_15;     
T_D1_15       = 170 - 3.33 * t_15;

t_step_25     = 0.0125;         
t_25          = 0:t_step_25:33;   
n_tsteps_25   = numel(t_25);     
n_sections_25 = n_tsteps_25;     
T_D1_25       = 170 - 3.33 * t_25;

t_step_3     = 0.015;         
t_3          = 0:t_step_3:33;   
n_tsteps_3   = numel(t_3);     
n_sections_3 = n_tsteps_3;     
T_D1_3       = 170 - 3.33 * t_3;

t_step_6     = 0.03;         
t_6          = 0:t_step_6:33;   
n_tsteps_6   = numel(t_6);     
n_sections_6 = n_tsteps_6;     
T_D1_6       = 170 - 3.33 * t_6;

t_step_8     = 0.04;         
t_8          = 0:t_step_8:33;   
n_tsteps_8   = numel(t_8);     
n_sections_8 = n_tsteps_8;     
T_D1_8       = 170 - 3.33 * t_8; 

%% Pre-allocate temperature and mask matrices
% Columns for D1 -- to end, and rows for time.
T_15 = zeros(n_tsteps_15);                    
mask_15 = zeros(n_tsteps_15, 'logical');      
for i = 1:667
	last_section_step_15 = n_sections_15 -i +1;
	T_15(i:end, i) = T_D1_15(1:last_section_step_15);
	mask_15(i:end, i) = true;
end

T_25 = zeros(n_tsteps_25);                    
mask_25 = zeros(n_tsteps_25, 'logical');      
for i = 1:400
	last_section_step_25 = n_sections_25 -i +1;
	T_25(i:end, i) = T_D1_25(1:last_section_step_25);
	mask_25(i:end, i) = true;
end

T_3 = zeros(n_tsteps_3);                    
mask_3 = zeros(n_tsteps_3, 'logical');      
for i = 1:333
	last_section_step_3 = n_sections_3 -i +1;
	T_3(i:end, i) = T_D1_3(1:last_section_step_3);
	mask_3(i:end, i) = true;
end

T_6 = zeros(n_tsteps_6);                    
mask_6 = zeros(n_tsteps_6, 'logical');      
for i = 1:167
	last_section_step_6 = n_sections_6 -i +1;
	T_6(i:end, i) = T_D1_6(1:last_section_step_6);
	mask_6(i:end, i) = true;
end

T_8 = zeros(n_tsteps_8);                    
mask_8 = zeros(n_tsteps_8, 'logical');
for i = 1:125
	last_section_step_8 = n_sections_8 -i +1;
	T_8(i:end, i) = T_D1_8(1:last_section_step_8);
	mask_8(i:end, i) = true;
end

%% Growth rate
T_threshold_15 = 120;   
G_15 = zeros(n_tsteps_15); 
above_threshold_15 = T_15 >= T_threshold_15;

G_above_15 = Goi * exp(-U ./ (Rg * (T_15 - 30))) .* exp(-Kgi .* (T_15 + C2K + Tm) ./ (2 * (T_15 + C2K).^2 .* (Tm - T_15 - C2K)));
G_below_15 = Goii * exp(-U ./ (Rg * (T_15 - 30))) .* exp(-Kgii .* (T_15 + C2K + Tm) ./ (2 * (T_15 + C2K).^2 .* (Tm - T_15 - C2K)));

G_15(above_threshold_15) = G_above_15(above_threshold_15);
G_15(~above_threshold_15) = G_below_15(~above_threshold_15);
G_15 = G_15 / 60;
G_15(~mask_15) = NaN;

T_threshold_25 = 120;   
G_25 = zeros(n_tsteps_25); 
above_threshold_25 = T_25 >= T_threshold_25;

G_above_25 = Goi * exp(-U ./ (Rg * (T_25 - 30))) .* exp(-Kgi .* (T_25 + C2K + Tm) ./ (2 * (T_25 + C2K).^2 .* (Tm - T_25 - C2K)));
G_below_25 = Goii * exp(-U ./ (Rg * (T_25 - 30))) .* exp(-Kgii .* (T_25 + C2K + Tm) ./ (2 * (T_25 + C2K).^2 .* (Tm - T_25 - C2K)));

G_25(above_threshold_25) = G_above_25(above_threshold_25);
G_25(~above_threshold_25) = G_below_25(~above_threshold_25);
G_25 = G_25 / 60;
G_25(~mask_25) = NaN;

T_threshold_3 = 120;   
G_3 = zeros(n_tsteps_3); 
above_threshold_3 = T_3 >= T_threshold_3;

G_above_3 = Goi * exp(-U ./ (Rg * (T_3 - 30))) .* exp(-Kgi .* (T_3 + C2K + Tm) ./ (2 * (T_3 + C2K).^2 .* (Tm - T_3 - C2K)));
G_below_3 = Goii * exp(-U ./ (Rg * (T_3 - 30))) .* exp(-Kgii .* (T_3 + C2K + Tm) ./ (2 * (T_3 + C2K).^2 .* (Tm - T_3 - C2K)));

G_3(above_threshold_3) = G_above_3(above_threshold_3);
G_3(~above_threshold_3) = G_below_3(~above_threshold_3);
G_3 = G_3 / 60;
G_3(~mask_3) = NaN;

T_threshold_6 = 120;   
G_6 = zeros(n_tsteps_6); 
above_threshold_6 = T_6 >= T_threshold_6;

G_above_6 = Goi * exp(-U ./ (Rg * (T_6 - 30))) .* exp(-Kgi .* (T_6 + C2K + Tm) ./ (2 * (T_6 + C2K).^2 .* (Tm - T_6 - C2K)));
G_below_6 = Goii * exp(-U ./ (Rg * (T_6 - 30))) .* exp(-Kgii .* (T_6 + C2K + Tm) ./ (2 * (T_6 + C2K).^2 .* (Tm - T_6 - C2K)));

G_6(above_threshold_6) = G_above_6(above_threshold_6);
G_6(~above_threshold_6) = G_below_6(~above_threshold_6);
G_6 = G_6 / 60;
G_6(~mask_6) = NaN;

T_threshold_8 = 120;   
G_8 = zeros(n_tsteps_8); 
above_threshold_8 = T_8 >= T_threshold_8;

G_above_8 = Goi * exp(-U ./ (Rg * (T_8 - 30))) .* exp(-Kgi .* (T_8 + C2K + Tm) ./ (2 * (T_8 + C2K).^2 .* (Tm - T_8 - C2K)));
G_below_8 = Goii * exp(-U ./ (Rg * (T_8 - 30))) .* exp(-Kgii .* (T_8 + C2K + Tm) ./ (2 * (T_8 + C2K).^2 .* (Tm - T_8 - C2K)));

G_8(above_threshold_8) = G_above_8(above_threshold_8);
G_8(~above_threshold_8) = G_below_8(~above_threshold_8);
G_8 = G_8 / 60;
G_8(~mask_8) = NaN;

%% Crystal size
% delta_Rp: Change (increase) in radius of one crystal at each time step
% This is the same for all sizes of crystals (i.e. existing crystals, and
% those nucleated at the current time step).delta_Rp_15 = G_15 * t_step_15;
delta_Rp_25= G_25 * t_step_25;
delta_Rp_3 = G_3 * t_step_3;
delta_Rp_6 = G_6 * t_step_6;
delta_Rp_8 = G_8 * t_step_8;

%% This is the section where MATLAB runs out of memory
% Rp: Current radius of one particle (crystal)
% Nucleation time (tn) is the page number i.e. Rp(t, D, tn) is the radius
% of one particle in section D, at time t, nucleated at time tn.

Rp_15 = zeros(n_tsteps_15, n_sections_15, n_tsteps_15);
for tn = 1:n_tsteps_15
    Rp_15(tn:end, :, tn) = cumsum(delta_Rp_15(tn:end, :), 1, 'omitnan');
end

Rp_25 = zeros(n_tsteps_25, n_sections_25, n_tsteps_25);
for tn = 1:n_tsteps_25
    Rp_25(tn:end, :, tn) = cumsum(delta_Rp_25(tn:end, :), 1, 'omitnan');
end

Rp_3 = zeros(n_tsteps_3, n_sections_3, n_tsteps_3);
for tn = 1:n_tsteps_3
    Rp_3(tn:end, :, tn) = cumsum(delta_Rp_3(tn:end, :), 1, 'omitnan');
end

Rp_6 = zeros(n_tsteps_6, n_sections_6, n_tsteps_6);
for tn = 1:n_tsteps_6
    Rp_6(tn:end, :, tn) = cumsum(delta_Rp_6(tn:end, :), 1, 'omitnan');
end

Rp_8 = zeros(n_tsteps_8, n_sections_8, n_tsteps_8);
for tn = 1:n_tsteps_8
    Rp_8(tn:end, :, tn) = cumsum(delta_Rp_8(tn:end, :), 1, 'omitnan');
end

% Vp: Volume of one particle (crystal) at each time step
Vp_15 = 4/3 * pi * Rp_15.^3;
Vp_25 = 4/3 * pi * Rp_25.^3;
Vp_3 = 4/3 * pi * Rp_3.^3;
Vp_6 = 4/3 * pi * Rp_6.^3;
Vp_8 = 4/3 * pi * Rp_8.^3;

%% PLA nucleation
% Total number of particles: not number of nucleated particles at current
% time step
N_15 =  No_15 * exp(-(1/k) ./ ((T_15 + C2K) .* (Tm - T_15 - C2K)));
N_25 =  No_25 * exp(-(1/k) ./ ((T_25 + C2K) .* (Tm - T_25 - C2K)));
N_3 =  No_3 * exp(-(1/k) ./ ((T_3 + C2K) .* (Tm - T_3 - C2K)));
N_6 =  No_6 * exp(-(1/k) ./ ((T_6 + C2K) .* (Tm - T_6 - C2K)));
N_8 =  No_8 * exp(-(1/k) ./ ((T_8 + C2K) .* (Tm - T_8 - C2K)));

% Number should be zero where T = 0. 
% This is becase where T = 0, it should actually equal NaN instead,
% but this yields an error in the calculation of G.
% Hence use mask afterwards to remove incorrect values.
N_15(~mask_15) = 0;
N_25(~mask_25) = 0;
N_3(~mask_3) = 0;
N_6(~mask_6) = 0;
N_8(~mask_8) = 0;

%% Accumulated volume
V_15 = zeros(n_tsteps_15, n_sections_15, n_tsteps_15);
V_25 = zeros(n_tsteps_25, n_sections_25, n_tsteps_25);
V_3 = zeros(n_tsteps_3, n_sections_3, n_tsteps_3);
V_6 = zeros(n_tsteps_6, n_sections_6, n_tsteps_6);
V_8 = zeros(n_tsteps_8, n_sections_8, n_tsteps_8);

% V: Total volume of all crystals
% V(t, D, tn) is the volume of all crystals nucleated at time tn, in
% section D1, at time t.

for tn = 1:n_tsteps_15
    V_15(:, :, tn) = N_15(tn, :) .* Vp_15(:, :, tn);
end

for tn = 1:n_tsteps_25
    V_25(:, :, tn) = N_25(tn, :) .* Vp_25(:, :, tn);
end

for tn = 1:n_tsteps_3
    V_3(:, :, tn) = N_3(tn, :) .* Vp_3(:, :, tn);
end

for tn = 1:n_tsteps_6
    V_6(:, :, tn) = N_6(tn, :) .* Vp_6(:, :, tn);
end

for tn = 1:n_tsteps_8
    V_8(:, :, tn) = N_8(tn, :) .* Vp_8(:, :, tn);
end

%% Total volume accumulated since t = 0.
% Vtot: total volume of all particles
                                
Vtot_nt_15 = sum(V_15, 3);        
Vtot_sec_15 = sum(Vtot_nt_15, 2); 
Vtot_t_15 = sum(Vtot_nt_15, 1);   
Vtot_all_15 = sum(Vtot_t_15);

Vtot_nt_25 = sum(V_25, 3);        
Vtot_sec_25 = sum(Vtot_nt_25, 2); 
Vtot_t_25 = sum(Vtot_nt_25, 1);   
Vtot_all_25 = sum(Vtot_t_25);
                            
Vtot_nt_3 = sum(V_3, 3);        
Vtot_sec_3 = sum(Vtot_nt_3, 2); 
Vtot_t_3 = sum(Vtot_nt_3, 1);   
Vtot_all_3 = sum(Vtot_t_3);

Vtot_nt_6 = sum(V_6, 3);        
Vtot_sec_6 = sum(Vtot_nt_6, 2); 
Vtot_t_6 = sum(Vtot_nt_6, 1);   
Vtot_all_6 = sum(Vtot_t_6);

Vtot_nt_8 = sum(V_8, 3);        
Vtot_sec_8 = sum(Vtot_nt_8, 2); 
Vtot_t_8 = sum(Vtot_nt_8, 1);   
Vtot_all_8 = sum(Vtot_t_8);   
                           