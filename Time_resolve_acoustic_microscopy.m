clearvars; close all; clc;

%% Frequency and transducer parameters
f0 = 50e6;
focal_length = 2.1e-3; %% scalling the focal_length from 12.8 mm to 2.1 mm for fist test          
half_aperture_angle = 14.4*pi/180; %% using the same angle that presented in the paper and convert the degree to radians by pi/180
aperture_diameter = 2*focal_length*sin(half_aperture_angle); %% calculate the aparture diameter based on the new facal length and the angle (sin(theta) = radius diameter/facal lenght) 

%% Materials  properties
c_water = 1481;  rho_water = 1000;
c_steel = 5900;  rho_steel = 7900;

% Attenuation properties (dB/cm))
% alpha_power = 2;             % k-Wave requires one power for the whole grid
% a_water = 2.2e-4;            % Water coefficient
% a_steel = 2.6e-4;            % Steel coefficient (from paper after convert from Np/m to dB/cm * 8.686)

%% Geometry
water_thickness_up   = 2.1e-3;
%h_film               = 200e-6;
h_film               = 230e-6; % I changed it to 230 um t have exactly the same the paper
water_thickness_down = 1e-3;

%% Grid
dx = 5e-6; %% we choose 5e-6 to have good soving the resplution but still the 1e-6 will be better but has slower of the simulation time
dy = dx;

total_thickness = water_thickness_up + h_film + water_thickness_down;
Nx = ceil(total_thickness/dx);
%Ny = ceil(4*aperture_diameter/dy);
Ny = ceil(2*aperture_diameter/dy); % I reduced the width more to reduce the simulation time
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% Calculate film positions
film_front = round(water_thickness_up/dx);
film_back  = film_front + round(h_film/dx);
film_front_before = round(1.9e-3/dx); %% Start the scan from here! this to make the scan start before the the front surface
fprintf('Film front: grid %d (%.3f mm)\n', film_front, film_front*dx*1e3);
%% Transducer geometry in grid points
focalLength_gp   = round(focal_length/dx);
diameter_gp = round(aperture_diameter/dx);

%% makeArc requires odd diameter (grid points)
if mod(diameter_gp, 2) == 0
    diameter_gp = diameter_gp + 1;
end

%% Source pulse parameters
bw = 0.6;
tc = gauspuls('cutoff', f0, bw);

%% Scan settings 
%z_um_vec_scan =0:5:1805;
z_um_vec_scan =1200;              
Nz = numel(z_um_vec_scan);

%% Time reducing z-factor starts at 2 and decreases to 0.5 with changing the z_um_vec
z_time = linspace(2, 0.5, Nz);

%% Storage of 2D array (cell arrays because Nt changes with z)
rf_meas_lines = cell(Nz, 1);   %% measurement: water / thin steel / water
rf_ref_lines  = cell(Nz, 1);   %% reference: water / bulk steel
t_lines       = cell(Nz, 1);   %% time array for each scan

t_end_all = zeros(Nz,1);
arc_z_mm  = zeros(Nz,1);
focus_z_mm = zeros(Nz,1);

%% Loop over defocus values
for iz = 1:Nz

    % Current depth (meters)
    z_m = z_um_vec_scan(iz) * 1e-6;

    % This is the line you asked to vary
    start_scan_with_steps = film_front_before + round(z_m/dx);    %% start from zero to 1800 um

    % Calculate arc position 
    transducer_position = start_scan_with_steps - focalLength_gp;

    % Keep arc inside the grid 
    if transducer_position < 2
        transducer_position = 2;
    end

    % Transducer position and focus position
    arc_pos   = [transducer_position, Ny/2];
    focus_pos = [arc_pos(1) + focalLength_gp, Ny/2];

    % For printing the result position in mm with each loop 
    arc_z_mm(iz)   = arc_pos(1) * dx * 1e3;
    focus_z_mm(iz) = focus_pos(1) * dx * 1e3;

    % Expected echo times (useful for choosing safe t_end)
    dist_to_front = (film_front - arc_pos(1)) * dx;          %% meters
    t_s1_exp = 2*dist_to_front / c_water;                    %% seconds (front echo)
    t_s2_exp = t_s1_exp + 2*h_film/c_steel;                  %% seconds (back echo)

    % Time array (your requested rule + safety so echoes fit)
    t_end_user = z_time(iz) * water_thickness_up/c_water + 2*h_film/c_steel; 
    t_end = max(t_end_user, t_s2_exp + 0.25e-6);             %% ensure S2 is inside if the t_end always bigger then t_s2_echo we can terminate the t_s2_echo
    t_end_all(iz) = t_end;

    kgrid.makeTime(max([c_water, c_steel]), 0.3, t_end);
    dt = kgrid.dt;
    t  = kgrid.t_array;
    Nt = numel(t);

    % Create source time series with correct length Nt (Center the pulse at tc so it starts near t=0)
    src_sig_full = gauspuls(t - tc, f0, bw);

    % Create arc mask
    transducer_mask = makeArc([Nx, Ny], arc_pos, focalLength_gp, diameter_gp, focus_pos);

    % Source and sensor (same for both cases)
    source.p_mask = transducer_mask; % the mask
    source.p      = src_sig_full; % the type of the wave

    sensor.mask   = transducer_mask; % the mask
    sensor.record = {'p'}; % the type of the recording here pressure can be change to other availabe values

    % water / thin steel film / water
   
    medium_meas.sound_speed = c_water * ones(Nx, Ny);
    medium_meas.density     = rho_water * ones(Nx, Ny);
    %medium_meas.alpha_coeff = a_water * ones(Nx, Ny); % <--- Add for atttnuation
    %medium_meas.alpha_power = alpha_power;           % <--- Add for atenuation
    medium_meas.sound_speed(film_front:film_back, :) = c_steel;
    medium_meas.density(film_front:film_back, :)     = rho_steel;
    %medium_meas.alpha_coeff(film_front:film_back, :) = a_steel; % <--- Add for atenuation

    %input_args = {'PMLSize', 10, 'PMLInside', false, 'PlotSim', false, 'DataCast', 'gpuArray-single'};
    input_args = {'PMLSize', 10, 'PMLInside', false, 'PlotSim', true, 'DataCast', 'single'};
    sensor_data_meas = kspaceFirstOrder2D(kgrid, medium_meas, source, sensor, input_args{:});

    rf_meas = sum(sensor_data_meas.p, 1);

    % water / bulk steel 
    
    medium_ref.sound_speed = c_water * ones(Nx, Ny);
    medium_ref.density     = rho_water * ones(Nx, Ny);
    %medium_ref.alpha_coeff = a_water * ones(Nx, Ny); % <--- Add for atttnuation
    %medium_ref.alpha_power = alpha_power;           % <--- Add for atttnuation
    medium_ref.sound_speed(film_front:end, :) = c_steel;
    medium_ref.density(film_front:end, :)     = rho_steel;
    %medium_ref.alpha_coeff(film_front:end, :) = a_steel; % <--- Add for atttnuation

    sensor_data_ref = kspaceFirstOrder2D(kgrid, medium_ref, source, sensor, input_args{:});

    rf_ref = sum(sensor_data_ref.p, 1);

    % Store
    rf_meas_lines{iz} = rf_meas;
    rf_ref_lines{iz}  = rf_ref;
    t_lines{iz}       = t;

    fprintf('Scan %d/%d | z=%4d um | arc=%.3f mm | focus=%.3f mm | Nt=%d | t_end=%.3f us\n', ...
        iz, Nz, z_um_vec_scan(iz), arc_z_mm(iz), focus_z_mm(iz), Nt, t_end*1e6);

end
%% -------------------- process and filter the simulation data -------------------------------

%% Build V(z,t) as 2D arrays (pad to Nt_max)

Nt_all = cellfun(@numel, t_lines);
Nt_max = max(Nt_all);

% Use dt from first scan to create a common time axis (seconds)
dt0 = t_lines{1}(2) - t_lines{1}(1);
t_common = (0:Nt_max-1) * dt0;

Vmeas = zeros(Nz, Nt_max);
Vref  = zeros(Nz, Nt_max);

for iz = 1:Nz
    rfM = rf_meas_lines{iz};
    rfR = rf_ref_lines{iz};
    Nt_i = numel(rfM);

    Vmeas(iz, 1:Nt_i) = rfM;
    %Vref(iz, 1:Nt_i)  = rfR;
end

%% filter the sinal from the excitation pulse
t_cut = 0.5e-6;   

for iz = 1:Nz
    rfM = rf_meas_lines{iz};
    %rfR = rf_ref_lines{iz};
    t_i = t_lines{iz};

    Nt_i = numel(t_i);

    % gate: keep only samples at/after 0.3 us
    gate = (t_i >= t_cut);            % logical vector length Nt_i

    rfM_g = rfM .* gate;
    %rfR_g = rfR .* gate;

    Vmeas(iz, 1:Nt_i) = rfM_g;
    %Vref(iz,  1:Nt_i) = rfR_g;
end
%% plot for testing the results:
iz = 1;  % choose which z index to plot (1..Nz)

figure('Color','w'); 
plot(t_common*1e6, Vmeas(iz,:), 'LineWidth', 1.5);
grid on;
xlabel('Time [\mus]');
ylabel('V_{meas} [a.u.]');
title(sprintf('V_{meas}(t) at z = %d \\mum (iz=%d)', z_um_vec_scan(iz), iz));
%% select the maxximum and make the window that start before and after 0.5 us

win_half = 0.3e-6;                 % 0.5 us before/after peak
dt0 = t_common(2) - t_common(1);   % common dt used in V(z,t)

Vmeas_win = zeros(size(Vmeas));
Vref_win  = zeros(size(Vref));

% (optional) store peak indices/times for debugging
ipk_meas = zeros(Nz,1);  tpk_meas = zeros(Nz,1);
ipk_ref  = zeros(Nz,1);  tpk_ref  = zeros(Nz,1);

for iz = 1:Nz

    
    vM = Vmeas(iz,:);
    [~, iM] = max(abs(vM));
    ipk_meas(iz) = iM;
    tpk_meas(iz) = t_common(iM);

    i1M = max(1, iM - round(win_half/dt0));
    i2M = min(numel(t_common), iM + round(win_half/dt0));
    Vmeas_win(iz, i1M:i2M) = vM(i1M:i2M);

    
    vR = Vref(iz,:);
    [~, iR] = max(abs(vR));
    ipk_ref(iz) = iR;
    tpk_ref(iz) = t_common(iR);

    i1R = max(1, iR - round(win_half/dt0));
    i2R = min(numel(t_common), iR + round(win_half/dt0));
    Vref_win(iz, i1R:i2R) = vR(i1R:i2R);

end

%% Plot first 4 traces from MEAS and REF in separate figures

nPlot = min(4, Nz);


figure('Name','MEAS: first 4 traces','Color','w');
for k = 1:nPlot
    subplot(2,2,k);
    idx = Nz - nPlot + k;
    plot(t_common*1e6, Vmeas_win(idx,:), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [\mus]');
    ylabel('V_{meas}');
    title(sprintf('MEAS trace %d (z = %d \\mum)', k, z_um_vec_scan(idx)));
end


figure('Name','REF: first 4 traces','Color','w');
for k = 1:nPlot
    subplot(2,2,k);
    idx = Nz - nPlot + k;
    plot(t_common*1e6, Vref_win(idx,:), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [\mus]');
    ylabel('V_{ref}');
    title(sprintf('REF trace %d (z = %d \\mum)', k, z_um_vec_scan(idx)));
end

%% -------------------- Spectrum analysis -------------------------------
% 
% %% --- Calculate Reflection Spectrum (R_measured) ---
% % Material impedances
% Zw = rho_water * c_water;
% Zs = rho_steel * c_steel;
% R_bulk_ref = (Zs - Zw) / (Zs + Zw); % Bulk reflection coefficient for normalization
% 
% % Axes
% fs = 1/dt0;
% Nt_win = size(Vmeas_win, 2);
% freq_vec = (0:Nt_win-1) * (fs / Nt_win);
% z_m_vec = z_um_vec_scan * 1e-6; 
% 
% % Temporal FFT
% SpecMeas = fft(Vmeas_win, [], 2); 
% SpecRef  = fft(Vref_win, [], 2);
% 
% % Extract frequency indices for 30-60 MHz
% fit_idx = find(freq_vec >= 28e6 & freq_vec <= 60e6);
% f_fit = freq_vec(fit_idx);
% 
% % Calculate R_measured using phase-corrected summation (Eq. 6)
% R_measured_full = zeros(1, Nt_win);
% for i = 1:length(fit_idx)
%     idx = fit_idx(i);
%     f = freq_vec(idx);
%     k0 = 2 * pi * f / c_water; % Wavenumber in water
% 
%     % Phase correction for transducer motion
%     phase_corr = exp(-2j * k0 * z_m_vec(:)); 
% 
%     % Spatial integration (Summation over z)
%     num = sum(SpecMeas(:, idx) .* phase_corr);
%     den = sum(SpecRef(:, idx) .* phase_corr);
% 
%     R_measured_full(idx) = (num / den) * R_bulk_ref;
% end
% 
% % Now R_measured exists for the fitting code
% R_measured = R_measured_full; 
% R_exp_mag = abs(R_measured(fit_idx));
% %% --- Inverse Algorithm (Fitting) ---
% 
% % Initial Guesses (based on steel properties)
% % x = [ZN (impedance ratio), h_bar (normalized TOF), alpha_bar (normalized atten)]
% % Note: h_bar = h_film / c_steel * 1e6 (assuming 1MHz ref)
% h_bar_guess = h_film * 1e6 / c_steel; 
% x0 = [Zs/Zw, h_bar_guess, 0.001]; 
% 
% % Optimization
% options = optimoptions('lsqnonlin', 'Display', 'final', 'FunctionTolerance', 1e-12);
% lb = [1, 0, 0];       % Lower bounds
% ub = [100, 1, 0.1];   % Upper bounds
% 
% [x_res, ~] = lsqnonlin(@(x) theory_error(x, f_fit, R_exp_mag), x0, lb, ub, options);
% 
% % --- Final Parameter Extraction ---
% ZN_final = x_res(1);
% h_bar_final = x_res(2);
% alpha_bar_final = x_res(3);
% 
% % 1. Velocity (c2)
% c_ext = (h_film * 1e6) / h_bar_final;
% 
% % 2. Density (rho2)
% rho_ext = (ZN_final * Zw) / c_ext;
% 
% % 3. Attenuation (alpha) at 50MHz
% f0_ref = 50e6;
% alpha_ext = alpha_bar_final * (2 * pi * f0_ref) / c_ext; 
% 
% fprintf('\n--- Extracted Results ---\n');
% fprintf('Velocity:  %.2f m/s (Actual: %.0f)\n', c_ext, c_steel);
% fprintf('Density:   %.2f kg/m^3 (Actual: %.0f)\n', rho_ext, rho_steel);
% fprintf('Attenuation at 50MHz: %.4f Np/m\n', alpha_ext);
% 
% %% --- Plot Comparison ---
% figure('Color','w');
% plot(f_fit*1e-6, R_exp_mag, 'r-.', 'MarkerSize', 12); hold on;
% % Calculate theoretical curve from results
% R_fit_mag = abs(calc_R_theory(x_res, f_fit));
% plot(f_fit*1e-6, R_fit_mag, 'b-', 'LineWidth', 2);
% xlabel('Frequency [MHz]','FontSize', 14, 'FontWeight', 'bold'); ylabel('Reflectiviy coeffecient','FontSize', 14, 'FontWeight', 'bold');
% set(gca, 'FontSize', 12);
% legend('Measured (Sim)', 'Theoretical Fit','FontSize', 12);
% ylim([0.3 1.2])
% 
% %% save all the results
% % out_name_1 = 'sim_results_all_oneUmStep_after_resize_1.mat';
% % save(out_name_1, '-v7.3');   
% 
% % Save simulation results
% % out_name = 'kwave_Vzt_meas_ref_part_oneUmStep.mat';
% % meta = struct();
% % meta.created_on = datestr(now);
% % meta.note = 'k-Wave 2D defocus scan: meas (water/thin steel/water) + ref (water/bulk steel)';
% % 
% % save(out_name, ...
% %     'meta', ...
% %     'rf_meas_lines', 'rf_ref_lines', 't_lines', 't_end_all', ...
% %     'Vmeas', 'Vref', 'Vmeas_win','Vref_win', 't_common', 'dt0', ...
% %     'f0', 'bw', 'tc', ...
% %     'dx', 'dy', 'Nx', 'Ny', ...
% %     'water_thickness_up', 'water_thickness_down', 'h_film', ...
% %     'c_water', 'rho_water', 'c_steel', 'rho_steel', ...
% %     'focal_length', 'half_aperture_angle', 'aperture_diameter', ...
% %     'film_front', 'film_back', 'film_front_before', ...
% %     'radius_gp', 'diameter_gp', ...
% %     'arc_z_mm', 'focus_z_mm', ...
% %     '-v7.3');
% % 
% % fprintf('Saved file: %s\n', out_name);
% %% --- Helper Functions ---
% function diff = theory_error(x, f, R_exp)
%     R_th = abs(calc_R_theory(x, f));
%     diff = R_th - R_exp;
% end
% 
% function R = calc_R_theory(x, f)
%     ZN = x(1);
%     h_bar = x(2);
%     alpha_bar = x(3);
% 
%     r12 = (ZN - 1) / (ZN + 1);
%     % phase term including normalized attenuation
%     phi = 2j * (2 * pi * f / 1e6) * h_bar * (1 + 1j * alpha_bar);
% 
%     % Three-layer reflection model (Eq. 2 in paper)
%     R = (r12 * (1 - exp(phi))) ./ (1 - r12^2 * exp(phi));
% end
% 
% %% ------------------- Reference (Truth) Values ---
% 
% % % Replace these with the actual values you used in your k-Wave medium properties
% % h_ref = 200e-6;      % Physical thickness of layer [m]
% % v_ref = 5650;          % Sound speed of layer (medium.sound_speed) [m/s]
% % rho_ref = 7900;        % Density of layer (medium.density) [kg/m^3]
% % alpha_ref = 7.528;     % Attenuation at 50 MHz [Np/m]
% % 
% % % --- Measured (Extracted) Values from your Inverse Algorithm ---
% % % These are the results obtained from your fitting code
% % h_meas = 200e-6;     % Extracted thickness 
% % v_meas = 6049;         % Extracted velocity
% % rho_meas = 6046;       % Extracted density
% % alpha_meas = 7.576;    % Extracted attenuation
% % 
% % % --- Calculate Relative Error Percentage ---
% % % Formula: |(Measured - Reference) / Reference| * 100
% % err_h = abs((h_meas - h_ref) / h_ref) * 100;
% % err_v = abs((v_meas - v_ref) / v_ref) * 100;
% % err_rho = abs((rho_meas - rho_ref) / rho_ref) * 100;
% % err_alpha = abs((alpha_meas - alpha_ref) / alpha_ref) * 100;
% % 
% % % --- Display Results ---
% % fprintf('--- Relative Error Results ---\n');
% % fprintf('Thickness Error: %.2f%%\n', err_h);
% % fprintf('Velocity Error:  %.2f%%\n', err_v);
% % fprintf('Density Error:   %.2f%%\n', err_rho);
% % fprintf('Attenuation Error: %.2f%%\n', err_alpha);
% 
% %% ------------------- several test that related to the 2D FFt using the sum method and second FFT method
% 
% 
% % %% --- Calculate Reflection Coefficient (Simplified) --- using summation the tomperal  FFt salution
% % fs = 1/dt0;
% % [Nz, Nt] = size(Vmeas_win);
% % freq_vec = (0:Nt-1) * (fs / Nt); % Full frequency vector in Hz
% % z_m_vec = z_um_vec_scan(:) * 1e-6;    % Ensure it's a column vector
% % 
% % % 1. Constants
% % Zw = rho_water * c_water;
% % Zs = rho_steel * c_steel;
% % R_bulk_ref = (Zs - Zw) / (Zs + Zw);
% % 
% % % 2. Temporal FFT (Full Range)
% % SpecMeas = fft(Vmeas_win, [], 2); 
% % SpecRef  = fft(Vref_win, [], 2);
% % 
% % % 3. Vectorized Phase Correction (The "Magic" Step)
% % % We create a matrix of wavenumbers for all frequencies
% % k0 = 2 * pi * freq_vec / c_water; 
% % 
% % % Matrix multiplication calculates the sum(Spec * exp(-2j*k0*z)) 
% % % for every frequency simultaneously.
% % phase_matrix = exp(-2j * z_m_vec * k0); % Size: [Nz x Nt]
% % num = sum(SpecMeas .* phase_matrix, 1);  % Sum across z-positions
% % den = sum(SpecRef .* phase_matrix, 1);   % Sum across z-positions
% % 
% % % 4. Final Ratio
% % R_meas_normal = (num ./ den) * R_bulk_ref;
% % 
% % %% --- Plotting with xlim ---
% % figure('Color', 'w');
% % plot(freq_vec * 1e-6, abs(R_meas_normal), 'LineWidth', 1.5);
% % grid on;
% % xlabel('Frequency [MHz]');
% % ylabel('|R_{meas}|');
% % title('Reflection Coefficient Spectrum');
% % 
% % % Focus only on the 30-60 MHz range
% % xlim([26 60]);
% % %% --- Calculate Reflection Coefficient using Spatial FFT --- Spatial FFT repated for comparison with sum method
% % 
% % % 1. Setup Parameters
% % Zw = rho_water * c_water;
% % Zs = rho_steel * c_steel;
% % R_bulk_ref = (Zs - Zw) / (Zs + Zw); % Bulk reference reflection
% % 
% % z_m = z_um_vec_scan * 1e-6;      % z-positions in meters
% % dz = z_m(2) - z_m(1);       % spatial step
% % Nz = size(Vmeas_win, 1);    % Number of z positions
% % Nt = size(Vmeas_win, 2);    % Number of time samples
% % fs = 1/dt0;                 % Sampling frequency
% % 
% % % Frequency axis
% % freq_vec = (0:Nt-1) * (fs / Nt);
% % f_idx = find(freq_vec >= 30e6 & freq_vec <= 60e6); % 30-60 MHz
% % 
% % % Initialize results
% % R_fft_method = zeros(size(f_idx));
% % 
% % % 2. Process each frequency
% % for i = 1:length(f_idx)
% %     idx = f_idx(i);
% %     f = freq_vec(idx);
% %     k0 = 2 * pi * f / c_water; % Wavenumber in water
% % 
% %     % Temporal FFT for this frequency (already a vector over z)
% %     Vz_meas_f = fft(Vmeas_win, [], 2); 
% %     Vz_ref_f  = fft(Vref_win, [], 2);
% % 
% %     % Extract the specific frequency component across all z
% %     Vz_meas = Vz_meas_f(:, idx);
% %     Vz_ref  = Vz_ref_f(:, idx);
% % 
% %     % Apply the phase correction e^(-2ik0z)
% %     % This aligns the signals relative to the focal plane
% %     phase_corr = exp(-2j * k0 * z_m(:));
% %     Vz_meas_corr = Vz_meas .* phase_corr;
% %     Vz_ref_corr  = Vz_ref .* phase_corr;
% % 
% %     % --- Spatial FFT Step ---
% %     % Transform from z-space to k-space
% %     % The k=0 component (normal incidence) is the first element of the FFT
% %     Fk_meas = fft(Vz_meas_corr); 
% %     Fk_ref  = fft(Vz_ref_corr);
% % 
% %     % Extract k=0 component
% %     num = Fk_meas(1); 
% %     den = Fk_ref(1);
% % 
% %     % Reflection coefficient
% %     R_fft_method(i) = (num / den) * R_bulk_ref;
% % end
% % 
% % %% --- Comparison Plot ---
% % % Note: To compare, the "Summation" method is: num = sum(Vz_meas_corr)
% % % Since fft(X)(1) == sum(X), the results will be identical.
% % 
% % figure('Color', 'w', 'Position', [100 100 800 400]);
% % plot(freq_vec(f_idx)*1e-6, abs(R_fft_method), 'r-', 'LineWidth', 2);
% % grid on;
% % xlabel('Frequency [MHz]');
% % ylabel('Reflectivity Magnitude |R|');
% % title('Reflectivity Spectrum using Spatial FFT (k=0)');
% % xlim([30 60]);
% % legend('Spatial FFT Method');