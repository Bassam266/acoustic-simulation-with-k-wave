%% Clean and close all
clearvars; close all; clc;

%%  Center frequency
f0 = 10e6;              

%% Acoustic properties 
c_water = 1480; rho_water = 1000; %taken from material list % c = (m.s-1); rho = kg.m-3
c_film   = 1450; rho_film   = 876; % oil properties taken from the paper % c = (m.s-1); rho = kg.m-3
%c_film   = 2600; rho_film   = 1200; % PMMA
c_glass = 5570; rho_glass = 2500; % glass teflon taken from material list % c = (m.s-1); rho = kg.m-3

%% Geometry
water_thickness_max = 1e-3; % (m)
glass_one_thickness = 1e-3; % (m)
glass_two_thickness = 1e-3; % (m)
film_thickness_list = [1.8e-6 5.3e-6 8.9e-6 12.4e-6 15.9e-6 19.5e-6]; % (m)
Nth = length(film_thickness_list);

%% Grid Calculation
dx = 1e-6;              % Steps resultion 
total_thickness = water_thickness_max + glass_one_thickness + glass_two_thickness + max(film_thickness_list);
Nx = 2^nextpow2(ceil(total_thickness / dx));
kgrid = kWaveGrid(Nx, dx);

% Time Array
t_end = (2*total_thickness / min([c_water, c_film, c_glass])) * 1.2;
kgrid.t_array = makeTime(kgrid, max([c_water, c_film, c_glass]), 0.3, t_end);
dt = kgrid.dt;


% Source (Gaussian Pulse)
source_position = 1;    % Source location index
source.p_mask = zeros(Nx, 1);
source.p_mask(source_position) = 1;
tp = -2/f0:dt:2/f0;
bw = 1.2; % bw = f_higher - f_lower/ f_center:: f_higher = 17 MHz, f_lower 4 MHz, f_center = ~10 MHz
input_signal = gauspuls(tp, f0, bw);
source.p = input_signal;

% Sensor
sensor.mask = source.p_mask;
sensor.record = {'p'};

%% reference signal simulation (back glass interface)
medium_ref.sound_speed = c_water * ones(Nx, 1);
medium_ref.density     = rho_water * ones(Nx, 1);

% Define only the first glass layer to get the back glass interface
% reflection. if i put air after it the simulation dose not work, no idea
% yet
g1_start = round(water_thickness_max / dx);
g1_end   = g1_start + round(glass_one_thickness / dx);
medium_ref.sound_speed(g1_start:g1_end) = c_glass;
medium_ref.density(g1_start:g1_end)     = rho_glass;

%input_args = {'PMLInside', false, 'PMLSize', 10, 'PlotSim', true, 'DataCast', 'single'};
input_args = {'PMLInside', false, 'PMLSize', 10, 'PlotSim', false, 'PlotFreq',300,'PlotLayout', true};
sensor_data_ref = kspaceFirstOrder1D(kgrid, medium_ref, source, sensor, input_args{:});
ref_signal = sensor_data_ref.p;

%% The glass–oil interface simulation

for i = 1:Nth
    current_oil = film_thickness_list(i);
    
    medium.sound_speed = c_water * ones(Nx, 1);
    medium.density     = rho_water * ones(Nx, 1);
    
    % Layer Boundaries
    g1_s = round(water_thickness_max / dx);
    g1_e = g1_s + round(glass_one_thickness / dx);
    oil_s = g1_e + 1;
    oil_e = oil_s + round(current_oil / dx);
    g2_s = oil_e + 1;
    g2_e = g2_s + round(glass_two_thickness / dx);
    
    medium.sound_speed(g1_s:g1_e) = c_glass; medium.density(g1_s:g1_e) = rho_glass;
    medium.sound_speed(oil_s:oil_e) = c_film; medium.density(oil_s:oil_e) = rho_film;
    medium.sound_speed(g2_s:g2_e) = c_glass; medium.density(g2_s:g2_e) = rho_glass;
    
    sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
    all_echo_signal(:,i) = sensor_data.p;
end

%% Windwing the signal and FFT analysis

% we will have the same windwo here because we are targeting the same echo
% position.
t_start = 1.78e-6;
t_end   = 2.01e-6;
time_window = (kgrid.t_array >= t_start & kgrid.t_array <= t_end).';

% sampling frequency
fs = 1/dt;

% FFT helper handle (gating + FFT) nice function for short code the
% function available in the repo ...
doFFT = @(signal) calc_fft2(signal .* time_window, fs);

% Reference FFT 
[freq_vec, specRef] = doFFT(ref_signal.');

% To make the film has the same size as reference 
Nfft = length(specRef);
specEcho   = zeros(Nfft, Nth);
reflecSpec = zeros(Nfft, Nth);

% Acoustic impedance and reflection coefficient
Zglass = rho_glass * c_glass;
Zwater = rho_water * c_water;

Rwg = (Zglass - Zwater) / (Zglass + Zwater);

% Loop over film echoes
for k = 1:Nth
    [~, specEcho(:,k)] = doFFT(all_echo_signal(:,k));
    
    % Reflection coefficient spectrum
    reflecSpec(:,k) = abs(specEcho(:,k) ./ specRef) * Rwg;
end

%% Ploting
figure('Name', 'Reflectivity Spectra Analysis');

% Time domain ref and film signal
subplot(2,1,1);
plot(kgrid.t_array*1e6, all_echo_signal(:,end), 'k'); hold on;
plot(kgrid.t_array*1e6, ref_signal, 'b--');
xregion(t_start*1e6, t_end*1e6, 'FaceColor', 'g', 'DisplayName', 'ref Window');
xregion(t_start*1e6, t_end*1e6, 'FaceColor', 'r', 'DisplayName', 'film Window');
xlabel('time (\mus)'); ylabel('pressure');
legend('film echo substrate film interface', 'Reference (Back first glass)');

% Reflectivity spectra
subplot(2,1,2);
for k = 1:Nth
    plot(freq_vec/1e6, reflecSpec(:,k), 'DisplayName', [num2str(film_thickness_list(k)*1e6) ' \mum']);
    hold on;
end
xlabel('frequency (MHz)'); ylabel('Reflectivity coefficients ');
xlim([6, 100]); ylim([0 1]);
grid on; legend('Location', 'bestoutside');

%% Seperate figure of the reflectivity
figure;
for k = 1:Nth
    plot(freq_vec/1e6, reflecSpec(:,k), 'DisplayName', [num2str(film_thickness_list(k)*1e6) ' \mum'], 'LineWidth', 2);
    hold on;
end
xlabel('Frequency (MHz)'); ylabel('Reflectivity coefficient');
xlim([6, 11]); ylim([0 1]);
grid on; legend('Location', 'bestoutside');
set(gca,"FontSize",14);
%% Thickness calcualtion
h_calc = zeros(length(freq_vec), Nth);
B = rho_film * (c_film^2);

% Ensure freq_vec is a column vector to match reflecSpec (Without that the freq_vec not the same!!)
if size(freq_vec, 1) == 1
    freq_vec = freq_vec'; 
end

for k = 1:Nth
    R = reflecSpec(:, k);
    term_R = sqrt( (R.^2) ./ (1 - R.^2) );
    
    % from 2.7 in Dwyer et al 2003
    % h = B / (pi * f * Z) * term_R
    denominator = pi .* freq_vec .* Zglass;
    h_calc(:, k) = (B ./ denominator) .* term_R;
end
%% Ploting
figure;

f_min = 1e6; 
f_max = f0;
f_idx = (freq_vec >= f_min & freq_vec <= f_max);

colors = lines(Nth); % Generate distinct colors

for k = 1:Nth
    % Plot the calculated thickness across the frequency range
    plot(freq_vec(f_idx)/1e6, h_calc(f_idx, k)*1e6, ...
        'Color', colors(k,:), 'LineWidth', 2, ...
        'DisplayName', ['expected thickness: ' num2str(film_thickness_list(k)*1e6) ' \mum']);
    
    % Calculate the average value 
    avg_h = mean(h_calc(f_idx, k)) * 1e6;
    fprintf('expected thickness: %.2f um | Calculated Avg: %.2f um\n', ...
        film_thickness_list(k)*1e6, avg_h);
    hold on;
end

% Formatting
xlabel('Frequency (MHz)');
ylabel('Thickness (\mum)');
xlim([6, 10]);
legend('Location', 'bestoutside');
set(gca, 'FontSize', 12);
%ylim([0, 5]);