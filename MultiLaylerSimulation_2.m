%% Ultrasound simulation - multilayer structure
% Water (4mm) -> Parylene-C or PMMA (10μm) -> Glass (1mm)
clear; close all; clc;
%% Parameters
f0 = 40e6;          

% Acoustic properties of materials
c_water = 1500;     
c_parylene = 2200;  
c_glass = 4200;     

rho_water = 1000;   
rho_parylene = 1280; 
rho_glass = 2500;   

% Layer geometry
water_thickness = .2e-3;     
parylene_thickness = 10e-6; 
glass_thickness = 1e-3;     
source_position = 1;       

%% Grid parameters
dx = 1e-6;

% Calculate minimum required grid size
total_thickness = water_thickness + parylene_thickness + glass_thickness;
min_grid_points = ceil(total_thickness / dx); %+ 100; % Add buffer % ceil is the function rounds each element of X to the nearest integer

% Ensure grid is large enough for all layers
Nx = max(min_grid_points);
kgrid = kWaveGrid(Nx, dx);

%% Medium properties 
medium.sound_speed = c_water * ones(Nx, 1);
medium.density = rho_water * ones(Nx, 1);

% Define layer boundaries
water_start = 1;
water_end = (round(water_thickness/dx));
parylene_start = water_end;
parylene_end = round(parylene_start + parylene_thickness/dx);
glass_start = parylene_end;
glass_end = (glass_start + round(glass_thickness / dx));


% Set Parylene-C layer properties

medium.sound_speed(parylene_start:parylene_end) = c_parylene;
medium.density(parylene_start:parylene_end) = rho_parylene;

% Set Glass layer properties

medium.sound_speed(glass_start:glass_end) = c_glass;
medium.density(glass_start:glass_end) = rho_glass;

% Attenuation setup (in dB/(MHz^y·cm))
medium.alpha_coeff = zeros(Nx,1);

medium.alpha_power = 1.5;
% Water
medium.alpha_coeff(water_start:water_end) = 0.0022;  % dB/(MHz^2·cm)

% Parylene-C
if parylene_start <= Nx && parylene_end <= Nx
    medium.alpha_coeff(parylene_start:parylene_end) = 0.2;  % dB/(MHz^1.5·cm)
end

% Glass
if glass_start <= Nx && glass_end <= Nx
    medium.alpha_coeff(glass_start:glass_end) = 0.1;  % dB/(MHz·cm)
end

%% Time array
% Calculate expected echo times for different interfaces
distance_water_parylene = (water_end+8) * dx; % - source_position
expected_echo_water_parylene = 2 * distance_water_parylene / c_water;

distance_parylene_water = (parylene_end - water_end) * dx;
expected_echo_parylene_water = 2 * distance_parylene_water / c_parylene;
expected_echo_parylene_glass = expected_echo_water_parylene + expected_echo_parylene_water;

distance_water_glass = (glass_end- source_position) * dx; %
expected_echo_water_glass = 2 * distance_water_glass / c_water;


% Use maximum sound speed for stability
c_max = max([c_water, c_parylene, c_glass]);
t_end = expected_echo_parylene_glass * 1.2; % Extended for multiple echoes
kgrid.t_array = makeTime(kgrid, c_max, 0.3, t_end);

%% Source definition
source.p_mask = zeros(Nx, 1);
source.p_mask(source_position) = 1;

tone_burst_cycles = 1;
input_signal = toneBurst(1/kgrid.dt, f0, tone_burst_cycles);
source.p = input_signal;

%% Sensor definition
sensor.mask = zeros(Nx, 1);
sensor.mask(source_position) = 1;
sensor.record = {'p'};

%% Run simulation with animation
fprintf('\nRunning simulation...\n');
input_args = {
    'PMLInside', false, ...
    'PMLSize', 10, ...
    'DataCast', 'single', ...
    'PlotSim', true, ...
    'PlotLayout', true ...
    };

sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
echo_signal_without_polymer = sensor_data.p;

fprintf('Simulation completed!\n');

%% run without polymer layer --> make polymer layer with water properties
 medium.sound_speed(parylene_start:parylene_end) = c_water;
 medium.density(parylene_start:parylene_end) = rho_water;
 medium.alpha_coeff(parylene_start:parylene_end) = 0.0022;  % dB/(MHz^2·cm)

fprintf('\nRunning simulation without polymer...\n');
input_args = {
    'PMLInside', false, ...
    'PMLSize', 10, ...
    'DataCast', 'single', ...
    'PlotSim', true, ...
    'PlotLayout', true ...
    };

sensor_data_nopoly = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
echo_signal_nopoly = sensor_data_nopoly.p;

fprintf('Simulation completed!\n');
%% plot difference

figure, nexttile, plot(kgrid.t_array,echo_signal_without_polymer),
hold on, plot(kgrid.t_array,echo_signal_nopoly), 
nexttile, plot(kgrid.t_array,echo_signal_without_polymer-echo_signal_nopoly)

%% downsample signals to match 500 MHz
downsample_factor = round(1 / (500e6 * kgrid.dt));
echo_signal_ds = downsample(echo_signal_without_polymer, downsample_factor);
echo_signal_nopoly_ds = downsample(echo_signal_nopoly, downsample_factor);
ds_time = downsample(kgrid.t_array, downsample_factor);

nexttile, plot(ds_time,echo_signal_ds),
hold on, plot(ds_time,echo_signal_nopoly_ds), 
nexttile, plot(ds_time,echo_signal_ds-echo_signal_nopoly_ds)
shg
return
%% edit plot
%% plot difference — styled 2x2 tiled layout
% % Consistent muted colors across both simulations
% colPoly   = [0.20, 0.40, 0.75];   % blue  -> with polymer
% colNoPoly = [0.85, 0.20, 0.20];   % red   -> without polymer
% colDiff   = [0.20, 0.20, 0.20];   % dark gray -> difference
% 
% figure;
% set(gcf, 'Color', 'w', 'Units', 'inches', 'Position', [1, 1, 10, 7]);
% tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% 
% % ---- Full-resolution simulation: overlaid signals ----
% nexttile;
% plot(kgrid.t_array * 1e6, echo_signal_without_polymer, ...
%     'Color', colPoly, 'LineWidth', 1.2); hold on;
% plot(kgrid.t_array * 1e6, echo_signal_nopoly, ...
%     'Color', colNoPoly, 'LineWidth', 1.2); hold off;
% xlabel('Time (\mus)',      'FontSize', 13, 'FontName', 'Arial');
% ylabel('Pressure (a.u.)',  'FontSize', 13, 'FontName', 'Arial');
% title ('Simulated echoes (full resolution)', 'FontSize', 13, 'FontName', 'Arial');
% lgd = legend('With polymer', 'Without polymer', 'Location', 'best');
% lgd.FontSize = 10; lgd.FontName = 'Arial'; lgd.Box = 'on';
% styleTimeAxes(gca);
% axis tight;
% 
% % ---- Full-resolution: difference ----
% nexttile;
% plot(kgrid.t_array * 1e6, echo_signal_without_polymer - echo_signal_nopoly, ...
%     'Color', colDiff, 'LineWidth', 1.2);
% xlabel('Time (\mus)',                'FontSize', 13, 'FontName', 'Arial');
% ylabel('Pressure difference (a.u.)', 'FontSize', 13, 'FontName', 'Arial');
% title ('Difference (with − without polymer)', 'FontSize', 13, 'FontName', 'Arial');
% styleTimeAxes(gca);
% axis tight;
% 
% % ---- Downsample to 500 MHz ----
% downsample_factor = round(1 / (500e6 * kgrid.dt));
% echo_signal_ds        = downsample(echo_signal_without_polymer, downsample_factor);
% echo_signal_nopoly_ds = downsample(echo_signal_nopoly,          downsample_factor);
% ds_time               = downsample(kgrid.t_array,               downsample_factor);
% 
% % ---- Downsampled: overlaid signals ----
% nexttile;
% plot(ds_time * 1e6, echo_signal_ds, ...
%     'Color', colPoly, 'LineWidth', 1.2); hold on;
% plot(ds_time * 1e6, echo_signal_nopoly_ds, ...
%     'Color', colNoPoly, 'LineWidth', 1.2); hold off;
% xlabel('Time (\mus)',      'FontSize', 13, 'FontName', 'Arial');
% ylabel('Pressure (a.u.)',  'FontSize', 13, 'FontName', 'Arial');
% title ('Simulated echoes (downsampled to 500 MHz)', 'FontSize', 13, 'FontName', 'Arial');
% lgd = legend('With polymer', 'Without polymer', 'Location', 'best');
% lgd.FontSize = 10; lgd.FontName = 'Arial'; lgd.Box = 'on';
% styleTimeAxes(gca);
% axis tight;
% 
% % ---- Downsampled: difference ----
% nexttile;
% plot(ds_time * 1e6, echo_signal_ds - echo_signal_nopoly_ds, ...
%     'Color', colDiff, 'LineWidth', 1.2);
% xlabel('Time (\mus)',                'FontSize', 13, 'FontName', 'Arial');
% ylabel('Pressure difference (a.u.)', 'FontSize', 13, 'FontName', 'Arial');
% title ('Difference (downsampled)', 'FontSize', 13, 'FontName', 'Arial');
% styleTimeAxes(gca);
% axis tight;
% 
% shg;
% 
% % ---- Local style helper ----
% function styleTimeAxes(ax)
%     ax.FontSize  = 11;
%     ax.FontName  = 'Arial';
%     ax.LineWidth = 1.2;
%     ax.Box       = 'on';
%     ax.TickDir   = 'in';
%     ax.XGrid     = 'on';
%     ax.YGrid     = 'on';
%     ax.GridAlpha = 0.2;
%     ax.GridColor = [0.5, 0.5, 0.5];
% end
