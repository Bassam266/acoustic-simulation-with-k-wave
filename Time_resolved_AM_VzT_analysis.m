%% Time-resolved acoustic microscopy: V(z,t) simulation + theoretical comparison
%
% Same physical idea as Time_resolve_acoustic_microscopy.m (defocus scan
% of a focused transducer over water/thin-film/water vs water/bulk
% reference), but this script actually runs the full pipeline end to end:
%
%   1) k-Wave defocus scan  ->  V(z,t) for "measured" (thin film) and
%      "reference" (bulk) cases
%   2) Temporal FFT + phase-corrected spatial summation over z (paper
%      Eq. 6)  ->  simulated reflectivity spectrum R_sim(f)
%   3) Fit the three-layer analytical reflection model (paper Eq. 2) to
%      R_sim(f) with lsqnonlin  ->  extracted thickness/velocity/
%      density/attenuation
%   4) Plot simulated spectrum vs. the fitted theoretical curve and
%      report the extraction error against the known (ground-truth)
%      medium properties used in the simulation.
%
% SCAN_MODE below switches between a coarse, few-point scan (fast, for
% testing the pipeline locally) and the fine step scan used for the
% actual paper-comparable result (run on a server/cluster).

clearvars; close all; clc;

%% ------------------------- Scan resolution switch -------------------------
% 'coarse' : few z-steps, quick local test of the whole pipeline
% 'fine'   : dense z-steps matching the paper's resolution, for the server
scan_mode = 'coarse';   % <-- change to 'fine' once verified locally

switch scan_mode
    case 'coarse'
        z_um_vec_scan = 0:1000:2000;   % ~10 points -> fast local test
    case 'fine'
        z_um_vec_scan = 0:5:1800;     % ~360 points -> paper-like resolution
    otherwise
        error('scan_mode must be ''coarse'' or ''fine''');
end
Nz = numel(z_um_vec_scan);

%% ------------------------- Frequency and transducer parameters -------------------------
f0 = 50e6;
focal_length = 2.1e-3;
half_aperture_angle = 14.4*pi/180;
aperture_diameter = 2*focal_length*sin(half_aperture_angle);

%% ------------------------- Materials properties -------------------------
c_water = 1481;  rho_water = 1000;
c_steel = 5900;  rho_steel = 7900;

%% ------------------------- Geometry -------------------------
water_thickness_up   = 2.1e-3;
h_film               = 230e-6;
water_thickness_down = 1e-3;

%% ------------------------- Grid -------------------------
dx = 5e-6;
dy = dx;

total_thickness = water_thickness_up + h_film + water_thickness_down;
Nx = ceil(total_thickness/dx);
Ny = ceil(2*aperture_diameter/dy);
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% ------------------------- Film position -------------------------
film_front = round(water_thickness_up/dx);
film_back  = film_front + round(h_film/dx);
film_front_before = round(1.9e-3/dx);   % scan starts before the front surface
fprintf('Film front: grid %d (%.3f mm)\n', film_front, film_front*dx*1e3);

%% ------------------------- Transducer geometry -------------------------
focalLength_gp = round(focal_length/dx);
diameter_gp    = round(aperture_diameter/dx);
if mod(diameter_gp, 2) == 0
    diameter_gp = diameter_gp + 1;
end

%% ------------------------- Source pulse -------------------------
bw = 0.6;
tc = gauspuls('cutoff', f0, bw);

%% ------------------------- Time-length shrink factor -------------------------
z_time = linspace(2, 0.5, Nz);

%% ------------------------- Storage -------------------------
rf_meas_lines = cell(Nz, 1);
rf_ref_lines  = cell(Nz, 1);
t_lines       = cell(Nz, 1);

t_end_all  = zeros(Nz,1);
arc_z_mm   = zeros(Nz,1);
focus_z_mm = zeros(Nz,1);

%% ------------------------- Loop over defocus values -------------------------
for iz = 1:Nz

    z_m = z_um_vec_scan(iz) * 1e-6;
    start_scan_with_steps = film_front_before + round(z_m/dx);

    transducer_position = start_scan_with_steps - focalLength_gp;
    if transducer_position < 2
        transducer_position = 2;
    end

    arc_pos   = [transducer_position, Ny/2];
    focus_pos = [arc_pos(1) + focalLength_gp, Ny/2];

    arc_z_mm(iz)   = arc_pos(1) * dx * 1e3;
    focus_z_mm(iz) = focus_pos(1) * dx * 1e3;

    dist_to_front = (film_front - arc_pos(1)) * dx;
    t_s1_exp = 2*dist_to_front / c_water;
    t_s2_exp = t_s1_exp + 2*h_film/c_steel;

    t_end_user = z_time(iz) * water_thickness_up/c_water + 2*h_film/c_steel;
    t_end = max(t_end_user, t_s2_exp + 0.25e-6);
    t_end_all(iz) = t_end;

    kgrid.makeTime(max([c_water, c_steel]), 0.3, t_end);
    dt = kgrid.dt;
    t  = kgrid.t_array;
    Nt = numel(t);

    src_sig_full = gauspuls(t - tc, f0, bw);

    transducer_mask = makeArc([Nx, Ny], arc_pos, focalLength_gp, diameter_gp, focus_pos);

    source.p_mask = transducer_mask;
    source.p      = src_sig_full;

    sensor.mask   = transducer_mask;
    sensor.record = {'p'};

    % --- measurement: water / thin steel film / water ---
    medium_meas.sound_speed = c_water * ones(Nx, Ny);
    medium_meas.density     = rho_water * ones(Nx, Ny);
    medium_meas.sound_speed(film_front:film_back, :) = c_steel;
    medium_meas.density(film_front:film_back, :)     = rho_steel;

    input_args = {'PMLSize', 10, 'PMLInside', false, 'PlotSim', false, 'DataCast', 'single'};
    sensor_data_meas = kspaceFirstOrder2D(kgrid, medium_meas, source, sensor, input_args{:});
    rf_meas = sum(sensor_data_meas.p, 1);

    % --- reference: water / bulk steel ---
    medium_ref.sound_speed = c_water * ones(Nx, Ny);
    medium_ref.density     = rho_water * ones(Nx, Ny);
    medium_ref.sound_speed(film_front:end, :) = c_steel;
    medium_ref.density(film_front:end, :)     = rho_steel;

    sensor_data_ref = kspaceFirstOrder2D(kgrid, medium_ref, source, sensor, input_args{:});
    rf_ref = sum(sensor_data_ref.p, 1);

    rf_meas_lines{iz} = rf_meas;
    rf_ref_lines{iz}  = rf_ref;
    t_lines{iz}       = t;

    fprintf('Scan %d/%d | z=%4d um | arc=%.3f mm | focus=%.3f mm | Nt=%d | t_end=%.3f us\n', ...
        iz, Nz, z_um_vec_scan(iz), arc_z_mm(iz), focus_z_mm(iz), Nt, t_end*1e6);
end

%% ============================ Post-processing ============================

%% --- Build V(z,t) as 2D arrays (pad to common Nt) ---
Nt_all = cellfun(@numel, t_lines);
Nt_max = max(Nt_all);

dt0 = t_lines{1}(2) - t_lines{1}(1);
t_common = (0:Nt_max-1) * dt0;

Vmeas = zeros(Nz, Nt_max);
Vref  = zeros(Nz, Nt_max);

t_cut = 0.5e-6;   % gate out the excitation pulse

for iz = 1:Nz
    rfM  = rf_meas_lines{iz};
    rfR  = rf_ref_lines{iz};
    t_i  = t_lines{iz};
    Nt_i = numel(t_i);

    gate = (t_i >= t_cut);

    Vmeas(iz, 1:Nt_i) = rfM .* gate;
    Vref(iz,  1:Nt_i) = rfR .* gate;
end

%% --- Window each trace around its peak echo ---
win_half = 0.3e-6;
dt0 = t_common(2) - t_common(1);

Vmeas_win = zeros(size(Vmeas));
Vref_win  = zeros(size(Vref));

for iz = 1:Nz
    vM = Vmeas(iz,:);
    [~, iM] = max(abs(vM));
    i1M = max(1, iM - round(win_half/dt0));
    i2M = min(numel(t_common), iM + round(win_half/dt0));
    Vmeas_win(iz, i1M:i2M) = vM(i1M:i2M);

    vR = Vref(iz,:);
    [~, iR] = max(abs(vR));
    i1R = max(1, iR - round(win_half/dt0));
    i2R = min(numel(t_common), iR + round(win_half/dt0));
    Vref_win(iz, i1R:i2R) = vR(i1R:i2R);
end

%% Quick sanity plots
nPlot = min(4, Nz);

figure('Name','MEAS: sample traces','Color','w');
for k = 1:nPlot
    subplot(2,2,k);
    idx = round(linspace(1, Nz, nPlot)); idx = idx(k);
    plot(t_common*1e6, Vmeas_win(idx,:), 'LineWidth', 1.5); grid on;
    xlabel('Time [\mus]'); ylabel('V_{meas}');
    title(sprintf('z = %d \\mum', z_um_vec_scan(idx)));
end

figure('Name','REF: sample traces','Color','w');
for k = 1:nPlot
    subplot(2,2,k);
    idx = round(linspace(1, Nz, nPlot)); idx = idx(k);
    plot(t_common*1e6, Vref_win(idx,:), 'LineWidth', 1.5); grid on;
    xlabel('Time [\mus]'); ylabel('V_{ref}');
    title(sprintf('z = %d \\mum', z_um_vec_scan(idx)));
end

%% ==================== Reflectivity spectrum (paper Eq. 6) ====================
Zw = rho_water * c_water;
Zs = rho_steel * c_steel;
R_bulk_ref = (Zs - Zw) / (Zs + Zw);   % bulk reflection coefficient (normalization)

fs = 1/dt0;
Nt_win = size(Vmeas_win, 2);
freq_vec = (0:Nt_win-1) * (fs / Nt_win);
z_m_vec = z_um_vec_scan(:) * 1e-6;

SpecMeas = fft(Vmeas_win, [], 2);
SpecRef  = fft(Vref_win,  [], 2);

fit_idx = find(freq_vec >= 28e6 & freq_vec <= 60e6);
f_fit = freq_vec(fit_idx);

k0_vec = 2*pi*f_fit/c_water;                       % [1 x Nf]
phase_matrix = exp(-2j * z_m_vec * k0_vec);        % [Nz x Nf]

num = sum(SpecMeas(:, fit_idx) .* phase_matrix, 1);
den = sum(SpecRef(:,  fit_idx) .* phase_matrix, 1);

R_measured = (num ./ den) * R_bulk_ref;
R_exp_mag  = abs(R_measured);

%% ==================== Fit theoretical model (paper Eq. 2) ====================
h_bar_guess = h_film * 1e6 / c_steel;
x0 = [Zs/Zw, h_bar_guess, 0.001];   % [ZN, h_bar, alpha_bar]

options = optimoptions('lsqnonlin', 'Display', 'final', 'FunctionTolerance', 1e-12);
lb = [1, 0, 0];
ub = [100, 1, 0.1];

[x_res, ~] = lsqnonlin(@(x) theory_error(x, f_fit, R_exp_mag), x0, lb, ub, options);

ZN_final     = x_res(1);
h_bar_final  = x_res(2);
alpha_bar_final = x_res(3);

c_ext   = (h_film * 1e6) / h_bar_final;
rho_ext = (ZN_final * Zw) / c_ext;
f0_ref  = 50e6;
alpha_ext = alpha_bar_final * (2*pi*f0_ref) / c_ext;

fprintf('\n--- Extracted Results (%s scan, Nz = %d) ---\n', scan_mode, Nz);
fprintf('Velocity:    %.2f m/s   (Actual: %.0f, error %.2f%%)\n', c_ext, c_steel, abs(c_ext-c_steel)/c_steel*100);
fprintf('Density:     %.2f kg/m^3 (Actual: %.0f, error %.2f%%)\n', rho_ext, rho_steel, abs(rho_ext-rho_steel)/rho_steel*100);
fprintf('Attenuation at 50MHz: %.4f Np/m\n', alpha_ext);
h_ext = h_bar_final * c_ext / 1e6;
fprintf('Thickness:   %.2f um    (Actual: %.2f um, error %.2f%%)\n', h_ext*1e6, h_film*1e6, abs(h_ext-h_film)/h_film*100);

%% --- Plot simulated spectrum vs. theoretical fit ---
figure('Color','w');
plot(f_fit*1e-6, R_exp_mag, 'r-.', 'LineWidth', 1.5); hold on;
R_fit_mag = abs(calc_R_theory(x_res, f_fit));
plot(f_fit*1e-6, R_fit_mag, 'b-', 'LineWidth', 2);
xlabel('Frequency [MHz]','FontSize', 14, 'FontWeight', 'bold');
ylabel('Reflectivity coefficient','FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
legend('Simulated (k-Wave)', 'Theoretical fit', 'FontSize', 12);
title(sprintf('Reflectivity spectrum: %s scan (Nz = %d, dz = %g \\mum)', ...
    scan_mode, Nz, mean(diff(z_um_vec_scan))));

%% --- Save results ---
out_name = sprintf('VzT_results_%s_Nz%d.mat', scan_mode, Nz);
save(out_name, 'z_um_vec_scan', 't_common', 'Vmeas', 'Vref', 'Vmeas_win', 'Vref_win', ...
    'freq_vec', 'f_fit', 'R_measured', 'x_res', 'c_ext', 'rho_ext', 'alpha_ext', 'h_ext', ...
    'f0', 'dx', 'dy', 'h_film', 'c_water', 'rho_water', 'c_steel', 'rho_steel', 'scan_mode', '-v7.3');
fprintf('Saved: %s\n', out_name);

%% ==================== Helper functions ====================
function diff = theory_error(x, f, R_exp)
    R_th = abs(calc_R_theory(x, f));
    diff = R_th - R_exp;
end

function R = calc_R_theory(x, f)
    ZN = x(1);
    h_bar = x(2);
    alpha_bar = x(3);

    r12 = (ZN - 1) / (ZN + 1);
    phi = 2j * (2*pi*f/1e6) * h_bar * (1 + 1j*alpha_bar);

    R = (r12 * (1 - exp(phi))) ./ (1 - r12^2 * exp(phi));
end
