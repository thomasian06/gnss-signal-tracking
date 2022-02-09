%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ian Thomas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

set(0, 'DefaultAxesLooseInset', [0,0,0,0])
set(0,'defaultAxesFontSize',14)

colors = get(gca, 'colororder');
close all;

addpath('../matlab');
addpath('data/');

rng(100);

%% Declaration of Parameters

Re = 6378.137e3; % semimajor axis of earth ellipsoid
f = 1/298.257223563; % flattening parameter of ellipsoid
e2 = 2*f - f^2; % eccentricity squared of ellipsoid
tol = 1e-10;
c = 299792458;

integration_periods = 1;

fs = 5e6;
ts = 1/fs;
fc = 1.023e6;
L1 = 1575.42e6; 
nc = 1023;
T = nc/fc * integration_periods;
nt = T * fs;
nsp = nt / integration_periods;
t = (0:(nt-1))' / fs;

fIF = 1.25e6;
fdmin = -5e3;
fdmax = 5e3;
fdspacing = 1/(2*T);
fd = fdmin:fdspacing:fdmax;

phi_t = @(t, fIF, fd) 2*pi*(fIF+fd).*t;

%% Simulated Data

prn = 31;
code = gen_prn_SV(prn);

simulated_CN0 = 60; % dB-Hz
simulated_SNR = 10^(simulated_CN0/10) / fs;
simulated_AWGN_sigma = sqrt(1/simulated_SNR);

simulated_doppler = (rand() * 2 - 1) * fdmax;
% simulated_doppler = 250;
simulated_carrier_phase = pi/2 * rand();
% simulated_carrier_phase = 0;
simulated_code_phase = 1023 * rand();
% simulated_code_phase = 500.5;

n_simulated_periods = 1000;
simulated_t = (0:(n_simulated_periods*nsp - 1)) / fs;

Sc = resample_digital(simulated_t, code, fc , simulated_code_phase);
Scarr = cos(phi_t(simulated_t, fIF, simulated_doppler)+...
    simulated_carrier_phase);

S = Sc .* Scarr + simulated_AWGN_sigma*randn(size(Sc));

gpsdata = S';


%% Acquisition

first = gpsdata(1:nt);
second = gpsdata(nt+1:2*nt);
third = gpsdata(2*nt+1:3*nt);

prns = 31;
SNR_acq = zeros(size(prns));
delays = zeros(size(prns));
dopplers = zeros(size(prns));
delay_inds = zeros(size(prns));
doppler_inds = zeros(size(prns));
noise_floors = zeros(size(prns));

% run coarse acquisition on all the prns
for i = 1:length(prns)
    code = gen_prn_SV(prns(i)) * 2 - 1;
    S_replica = resample_digital(t, code, fc, 0)';
    [lags, corrs] = cxcorr(first,...
        S_replica.*exp(1i*phi_t(t, fIF, fd)), nt);
    lags = lags(1:(nt/integration_periods));
    corrs = corrs(1:(nt/integration_periods), :);
    abscorrs = abs(corrs) / nt;
    [maxcorrs, maxidx] = max(abscorrs(:));
    maxcol = floor(maxidx / (nt/integration_periods)) + 1;
    maxrow = mod(maxidx - 1, (nt/integration_periods)) + 1;
    delays(i) = lags(maxrow);
    dopplers(i) = fd(maxcol);
    delay_inds(i) = maxrow;
    doppler_inds(i) = maxcol;
    SNR_acq(i) = 20*log10(max(abscorrs(:)) / mean(abscorrs(:)));
    noise_floors(i) = mean(abscorrs(:));
    
    figure(3*i-2)
    hold on;
    grid on;
    xlabel('Delay [Samples]')
    ylabel('Doppler [Hz]')
    zlabel('$S$')
    title(sprintf("PRN %d", prns(i)))
    [X, Y] = meshgrid(lags, fd);
    surface(X, Y, abscorrs'/max(abscorrs(:)), 'edgecolor', 'none')
    shading interp;
    view(35, 35);
    saveas(gcf, sprintf("figures/prn%d_acq_i%d_CN0%d", prns(i),...
        integration_periods, simulated_CN0), 'epsc')
    
    figure(3*i-1)
    hold on;
    grid on;
    xlabel('Delay [Samples]')
    ylabel('$S$')
    title(sprintf('$S$ vs. Delay, PRN %d', prns(i)));
    plot(lags, abscorrs(:, maxcol) / maxcorrs);
    saveas(gcf, sprintf("figures/prn%d_del_i%d_CN0%d", prns(i),...
        integration_periods, simulated_CN0), 'epsc')
    
    figure(3*i)
    hold on;
    grid on;
    xlabel('Doppler [Hz]')
    ylabel('$S$')
    title(sprintf('$S$ vs. Doppler, PRN %d', prns(i)));
    plot(fd, abscorrs(maxrow, :) / maxcorrs);
    saveas(gcf, sprintf("figures/prn%d_dop_i%d_CN0%d", prns(i),...
        integration_periods, simulated_CN0), 'epsc')
end

% use an acquisition SNR of 15 dB as the cutoff
cutoff = SNR_acq >= 15;
SNR_acq = SNR_acq(cutoff);
prns = prns(cutoff);
delays = delays(cutoff);
dopplers = dopplers(cutoff);
delay_inds = delay_inds(cutoff);
doppler_inds = doppler_inds(cutoff);
noise_floors = noise_floors(cutoff)';


%% Tracking

ns = length(gpsdata);
nb = floor(ns / nt);
epl_delay = 0.5; % chip difference between early, prompt, and late correlators

B_DLL = 0.2; % tracking loop bandwidths
B_PLL = 0.01/(2*pi);
B_FLL = 12;

gpsdata = gpsdata(1:(nb*nt));
gpsblocks = reshape(gpsdata, [nt, nb]);

prn_code_phases = zeros(length(prns), nb);
prn_carrier_phases = zeros(length(prns), nb);
prn_doppler_freqs = zeros(length(prns), nb);
prn_ecorrs = zeros(length(prns), nb);
prn_pcorrs = zeros(length(prns), nb);
prn_lcorrs = zeros(length(prns), nb);
prn_delta_thetas = zeros(length(prns), nb);
prn_delta_f = zeros(length(prns), nb);

for i = 1:length(prns)
    prn = prns(i);
    code_phase = nc - delays(i) / nsp * nc;
    carrier_phase = 0;
    doppler_freq = dopplers(i);
    
    code_phases = zeros(size(1:nb));
    carrier_phases = zeros(size(1:nb));
    doppler_freqs = zeros(size(1:nb));
    ecorrs = zeros(size(1:nb));
    pcorrs = zeros(size(1:nb));
    lcorrs = zeros(size(1:nb));
    delta_thetas = zeros(size(1:nb));
    delta_f = zeros(size(1:nb));
    
    IP1 = 0;
    QP1 = 0;
    
    IC_delta_theta = 0;
    
    for j = 1:nb
        
        block = gpsblocks(:, j);
        block = block - mean(block);
        
        code = gen_prn_SV(prn);
        
        carrier_replica = exp(-1i*(phi_t(t, fIF, doppler_freq)...
            + carrier_phase));
        
        % early
        ecode = resample_digital(t, code, fc + (1 + doppler_freq/L1),...
            -epl_delay + code_phase) * 2 - 1;
        Se = ecode' .* carrier_replica;
        ecorr = sum(block .* conj(Se)) / nt;
        ecorrs(j) = ecorr;
        
        % prompt
        pcode = resample_digital(t, code, fc + (1 + doppler_freq/L1),...
            code_phase) * 2 - 1;
        Sp = pcode' .* carrier_replica;
        pcorr = sum(block .* conj(Sp)) / nt;
        pcorrs(j) = pcorr;
        
        % late
        lcode = resample_digital(t, code, fc + (1 + doppler_freq/L1),...
            epl_delay + code_phase) * 2 - 1;
        Sl = lcode' .* carrier_replica;
        lcorr = sum(block .* conj(Sl)) / nt;
        lcorrs(j) = lcorr;
        
        % DLL
        dll_discriminator = -epl_delay*(abs(ecorr) - abs(lcorr)) / ...
            (abs(ecorr) + abs(lcorr));
        code_phase_error = B_DLL * dll_discriminator;
        code_phase = code_phase + code_phase_error;
        
        % PLL
        pll_discriminator = atan2(imag(pcorr), real(pcorr));
        delta_thetas(j) = pll_discriminator;
        IC_delta_theta = IC_delta_theta + pll_discriminator;
        carrier_phase_error = -B_PLL * IC_delta_theta;
        carrier_phase = carrier_phase + carrier_phase_error;

        % FLL
        IP2 = real(pcorr);
        QP2 = imag(pcorr);
        
        fll_cross = IP1.*QP2 - IP2.*QP1;
        fll_dot = IP1.*IP2 + QP1.*QP2;
        fll_discriminator = mean(atan2(fll_cross, fll_dot));
        
        IP1 = IP2;
        QP1 = QP2;
        
        doppler_freq = doppler_freq - B_FLL * fll_discriminator;
        doppler_freqs(j) = doppler_freq;
        delta_f(j) = fll_discriminator;
        
        f_code_adj = fc * (1+ doppler_freq/L1);
        code_phase = code_phase + f_code_adj*T;
        code_phases(j) = mod(code_phase, nc);
        carrier_phase = carrier_phase + 2*pi*(fIF + doppler_freq)*T;
        carrier_phases(j) = mod(carrier_phase, 2*pi);
        carrier_phases(j) = carrier_phase;
    end
    
    prn_code_phases(i, :) = code_phases;
    prn_carrier_phases(i, :) = carrier_phases;
    prn_doppler_freqs(i, :) = doppler_freqs;
    prn_ecorrs(i, :) = ecorrs;
    prn_pcorrs(i, :) = pcorrs;
    prn_lcorrs(i, :) = lcorrs;
    prn_delta_thetas(i, :) = delta_thetas;
    prn_delta_f(i, :) = delta_f;

end


%% Compare new code phase correlations to old ones

figure(10)
hold on;
grid on;
xlabel('Block')
ylabel('Correlation')
title('EPL Correlators vs. Block')
plot(1:nb, abs(pcorrs));
plot(1:nb, abs(ecorrs));
plot(1:nb, abs(lcorrs));
legend('Prompt', 'Early', 'Late', 'location', 'se')
saveas(gcf, sprintf("figures/epl_i%d_CN0%d", integration_periods,...
    simulated_CN0), 'epsc')

figure(11)
hold on;
grid on;
xlabel('Block')
ylabel('$\Delta\theta$')
title('Carrier Phase Misalignment vs. Block')
plot(1:nb, prn_delta_thetas, 'LineWidth', 1.5)
saveas(gcf, sprintf("figures/dtheta_i%d_CN0%d", integration_periods,...
    simulated_CN0), 'epsc')

figure(12)
hold on;
grid on;
xlabel('Block')
ylabel('Tracking SNR')
title('Tracking Autocorrelation SNR vs. Block')
SNR_trk = 20*log10(abs(prn_pcorrs)./noise_floors)';
plot(1:nb, SNR_trk, 'LineWidth', 1.5)
yline(SNR_acq, '--', 'LineWidth', 1.5, 'Color', colors(2, :))
legend('SNR', 'Acquisition SNR', 'location', 'se')
saveas(gcf, sprintf("figures/SNR_i%d_CN0%d", integration_periods,...
    simulated_CN0), 'epsc')

figure(14)
hold on;
grid on;
xlabel('Block')
ylabel('Doppler')
title('Tracked Doppler vs. Block')
plot(1:nb, prn_doppler_freqs, 'LineWidth', 1.5)
yline(simulated_doppler, '--r', 'LineWidth', 1.5)
legend('Tracked Doppler', 'True Doppler')
saveas(gcf, sprintf("figures/doppler_i%d_CN0%d", integration_periods,...
    simulated_CN0), 'epsc')


figure(15)
hold on;
grid on;
xlabel('Block')
ylabel('Code Phase')
title('Tracked Code Phase vs. Block')
plot(1:nb, prn_code_phases, 'LineWidth', 1.5)
yline(nc - delays(i) / nsp * nc, '--', 'LineWidth', 1.5,...
    'Color', colors(2, :))
yline(nc - (delays(i)+1) / nsp * nc, '--', 'LineWidth', 1.5,...
    'Color', colors(3, :))
yline(simulated_code_phase, '--r', 'LineWidth', 1.5)
legend('Code Phase', 'Sample Boundary', 'Sample Boundary',...
    'True Code Phase', 'location', 'se')
saveas(gcf, sprintf("figures/code_phase_i%d_CN0%d",...
    integration_periods, simulated_CN0), 'epsc')
