% DPD_sampling_analysis.m  --  DPD Basis Function Aliasing Analysis
% =============================================================================
% Question:
%   DPD polynomial order K 的 basis  y_k = x*|x|^(k-1) bandwidth = k*BW。
%   如果 Fs = M*BW 且 K > M，超出 Nyquist 的頻譜折回 → in-band aliasing。
%   這裡用 12X 當 golden reference，比較低 Fs 的 waveform error。
%
% Method:
%   1. Generate signal at Fs_ref = 720 MHz (36X, alias-free for all K ≤ 18)
%   2. Decimate to each Fs = [3,6,9,12]*BW
%   3. Compute y_k = x*|x|^(K-1) at each Fs
%   4. Upsample back to Fs_ref via interpft (ideal bandlimited interpolation)
%   5. Compare with golden (computed at Fs_ref): NMSE total + in-band
%
% PA model:
%   y = x * (a1 + a3*|x|^2 + a5*|x|^4 + a7*|x|^6)
%   Used for AM/AM + AM/PM characterization and ACLR context.
% =============================================================================

clear; clc; close all;

%% 1. Parameters
BW     = 20e6;
K_vec  = [1, 3, 5, 7, 9];
M_vec  = [3, 6, 9, 12];
nK     = length(K_vec);
nM     = length(M_vec);

Fs_ref = 720e6;                   % 36X, alias-free reference
N_ref  = 98304;                   % 32768*3, divisible by D = [12,6,4,3]
D_vec  = Fs_ref ./ (M_vec * BW); % decimation factors: [12, 6, 4, 3]

% PA coefficients (for characterization figure)
pa_c = [1.0, -0.3+0.1j, 0.08-0.05j, -0.02+0.01j];

% DPD regularization
lam = 1e-6;

fprintf('=== DPD Basis Function Aliasing Analysis ===\n');
fprintf('BW = %.0f MHz,  Fs_ref = %.0f MHz (%.0fX)\n', BW/1e6, Fs_ref/1e6, Fs_ref/BW);
fprintf('DPD orders K = %s\n', mat2str(K_vec));
fprintf('Oversampling M = %s  ->  Fs = %s MHz\n', ...
        mat2str(M_vec), mat2str(M_vec * BW/1e6));
fprintf('Decimation factors D = %s\n', mat2str(D_vec));
fprintf('N_ref = %d,  N at each Fs = %s\n\n', ...
        N_ref, mat2str(N_ref ./ D_vec));

%% 2. PA Characterization (Figure 1)
r_ax = linspace(0, 1.2, 1000)';
G_pa = pa_c(1) + pa_c(2)*r_ax.^2 + pa_c(3)*r_ax.^4 + pa_c(4)*r_ax.^6;

figure('Position', [80 80 1300 380]);
subplot(1,3,1);
plot(r_ax, abs(G_pa).*r_ax, 'b-', 'LineWidth', 2); hold on;
plot(r_ax, r_ax, 'k--', 'LineWidth', 1);
xlabel('|x_{in}|'); ylabel('|y_{out}|');
title('AM/AM'); grid on;
legend('PA', 'Linear', 'Location', 'northwest');

subplot(1,3,2);
plot(r_ax, angle(G_pa)*180/pi, 'r-', 'LineWidth', 2);
xlabel('|x_{in}|'); ylabel('Phase shift (deg)');
title('AM/PM'); grid on;

subplot(1,3,3);
gain_dB = 20*log10(abs(G_pa));
gain_dB(r_ax < 0.01) = NaN;
plot(r_ax, gain_dB, 'b-', 'LineWidth', 2); hold on;
yline(0, 'k--'); yline(-1, 'r:', 'P1dB');
xlabel('|x_{in}|'); ylabel('Gain (dB)');
title('Gain Compression'); grid on; ylim([-10 2]);

sgtitle('PA Model: y = x \cdot (a_1 + a_3|x|^2 + a_5|x|^4 + a_7|x|^6)');
saveas(gcf, 'DPD_sampling_PA_char.png');
fprintf('Saved DPD_sampling_PA_char.png\n');

%% 3. Signal Generation at Fs_ref
rng(42);
df    = Fs_ref / N_ref;
f_ref = ((-N_ref/2):(N_ref/2-1))' * df;
mask  = double(abs(f_ref) <= BW/2);
x_ref = ifft(ifftshift(mask .* exp(1j * 2*pi * rand(N_ref, 1))));
x_ref = x_ref / max(abs(x_ref));

rms_x  = sqrt(mean(abs(x_ref).^2));
fprintf('Signal: peak=1, RMS=%.4f, PAPR=%.2f dB\n\n', rms_x, 20*log10(1/rms_x));

%% 4. Main Analysis: Golden vs Aliased Basis Functions
NMSE_total  = NaN(nK, nM);
NMSE_inband = NaN(nK, nM);
psd_ref     = cell(nK, 1);       % golden PSD at Fs_ref
psd_up      = cell(nK, nM);      % upsampled aliased PSD

f_ref_MHz = f_ref / 1e6;
ib_mask   = abs(f_ref) <= BW/2;  % in-band mask at Fs_ref

fprintf('%-5s %-5s %-10s %12s %12s\n', ...
        'K', 'M', 'condition', 'NMSE_total', 'NMSE_inband');
fprintf('%s\n', repmat('-', 1, 46));

for ii = 1:nK
    K = K_vec(ii);

    % Golden: basis at Fs_ref (alias-free)
    y_ref = x_ref .* abs(x_ref).^(K-1);
    Y_ref = fftshift(fft(y_ref)) / N_ref;
    psd_ref{ii} = abs(Y_ref).^2;

    for jj = 1:nM
        M = M_vec(jj);
        D = D_vec(jj);

        % Decimate signal to Fs = M*BW
        x_lo = x_ref(1:D:end);

        % Compute basis at low Fs (aliased when K > M)
        y_lo = x_lo .* abs(x_lo).^(K-1);

        % Upsample to Fs_ref for comparison
        y_up = interpft(y_lo, N_ref);
        Y_up = fftshift(fft(y_up)) / N_ref;
        psd_up{ii, jj} = abs(Y_up).^2;

        % Total waveform NMSE
        err = y_up - y_ref;
        NMSE_total(ii, jj) = 10*log10(mean(abs(err).^2) / mean(abs(y_ref).^2));

        % In-band NMSE (within +/- BW/2)
        err_ib = sum(abs(Y_up(ib_mask) - Y_ref(ib_mask)).^2);
        ref_ib = sum(abs(Y_ref(ib_mask)).^2);
        if ref_ib < 1e-30
            NMSE_inband(ii, jj) = -Inf;
        else
            NMSE_inband(ii, jj) = 10*log10(err_ib / ref_ib);
        end

        % Tag
        if     K < M,  tag = 'OK';
        elseif K == M, tag = 'Edge';
        else,          tag = 'ALIASED';
        end

        fprintf('%-5d %-5d %-10s %12.2f %12.2f\n', ...
                K, M, tag, NMSE_total(ii,jj), NMSE_inband(ii,jj));
    end
end

%% 5. Figure 2 -- Spectral Grid (5 x 4)
%  Each subplot: golden spectrum (gray) vs aliased+reconstructed (blue)
figure('Position', [30 30 1600 1000]);
for ii = 1:nK
    K = K_vec(ii);
    ref_peak = max(psd_ref{ii});

    for jj = 1:nM
        M   = M_vec(jj);
        Fs  = M * BW;
        sp  = (ii - 1) * nM + jj;
        subplot(nK, nM, sp);

        % Golden (gray area)
        sdb_ref = 10*log10(psd_ref{ii} / ref_peak + 1e-30);
        area(f_ref_MHz, max(sdb_ref, -100), 'BaseValue', -100, ...
             'FaceColor', [0.85 0.85 0.85], 'EdgeColor', [0.6 0.6 0.6]);
        hold on;

        % Aliased+reconstructed (blue)
        sdb_up = 10*log10(psd_up{ii,jj} / ref_peak + 1e-30);
        plot(f_ref_MHz, sdb_up, 'b-', 'LineWidth', 0.8);

        % In-band edges
        xline(-BW/2/1e6, 'k--', 'LineWidth', 0.5);
        xline( BW/2/1e6, 'k--', 'LineWidth', 0.5);

        % Nyquist of this Fs (red dashed)
        xline(-Fs/2/1e6, 'r:', 'LineWidth', 0.8);
        xline( Fs/2/1e6, 'r:', 'LineWidth', 0.8);

        % Expected basis bandwidth (green dashed)
        xline(-K*BW/2/1e6, ':', 'Color', [0 0.6 0], 'LineWidth', 0.6);
        xline( K*BW/2/1e6, ':', 'Color', [0 0.6 0], 'LineWidth', 0.6);

        xlim([-150 150]); ylim([-100 5]); grid on;

        % Tag + NMSE annotation
        if K > M
            text(0.03, 0.88, sprintf('ALIASED\nNMSE_{ib}=%.1f dB', ...
                 NMSE_inband(ii,jj)), 'Units','normalized', ...
                 'Color',[0.8 0 0], 'FontWeight','bold', 'FontSize',6);
        elseif K == M
            text(0.03, 0.88, sprintf('EDGE\nNMSE_{ib}=%.1f dB', ...
                 NMSE_inband(ii,jj)), 'Units','normalized', ...
                 'Color',[0.8 0.5 0], 'FontWeight','bold', 'FontSize',6);
        else
            text(0.03, 0.88, sprintf('OK\nNMSE_{ib}=%.1f dB', ...
                 NMSE_inband(ii,jj)), 'Units','normalized', ...
                 'Color',[0 0.6 0], 'FontWeight','bold', 'FontSize',6);
        end

        if ii == 1
            title(sprintf('M=%d (Fs=%dMHz)', M, Fs/1e6), 'FontSize', 9);
        end
        if jj == 1
            ylabel(sprintf('K=%d\n(BW=%dMHz)', K, K*BW/1e6), ...
                   'FontSize', 8, 'FontWeight', 'bold');
        end
        if ii == nK, xlabel('MHz', 'FontSize', 7); end
        set(gca, 'FontSize', 6);
    end
end
sgtitle({'DPD Basis Aliasing: y_k = x \cdot |x|^{K-1}', ...
         'gray = golden (36X), blue = decimated+reconstructed, red: = Nyquist, green: = basis BW'});
saveas(gcf, 'DPD_sampling_spectra.png');
fprintf('\nSaved DPD_sampling_spectra.png\n');

%% 6. Figure 3 -- Total NMSE Heatmap (full bandwidth waveform error)
NMSE_disp = max(NMSE_total, -100);  % clamp for display

figure('Position', [100 100 640 500]);
imagesc(NMSE_disp);
colormap(parula); cb = colorbar;
ylabel(cb, 'NMSE_{total} (dB)');
caxis([-105, max(NMSE_disp(:)) + 3]);
set(gca, 'XTick', 1:nM, 'XTickLabel', ...
    arrayfun(@(m) sprintf('%dX', m), M_vec, 'UniformOutput', false));
set(gca, 'YTick', 1:nK, 'YTickLabel', ...
    arrayfun(@(k) sprintf('K=%d', k), K_vec, 'UniformOutput', false));
xlabel('Oversampling  M = Fs / BW');
ylabel('DPD Polynomial Order  K');
title('Total Waveform Error (dB): Basis at Fs vs Golden');

for ii = 1:nK
    for jj = 1:nM
        v = NMSE_total(ii, jj);
        if v < -100, s = '0'; else, s = sprintf('%.1f', v); end
        if v > -60, clr = 'w'; else, clr = 'k'; end
        text(jj, ii, s, 'HorizontalAlignment', 'center', ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', clr);
        if K_vec(ii) > M_vec(jj)
            rectangle('Position', [jj-0.5 ii-0.5 1 1], ...
                      'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
        end
    end
end
saveas(gcf, 'DPD_sampling_NMSE.png');
fprintf('Saved DPD_sampling_NMSE.png\n');

%% Figure 3b -- In-Band NMSE Heatmap
NMSE_ib_disp = max(NMSE_inband, -100);

figure('Position', [110 110 640 500]);
imagesc(NMSE_ib_disp);
colormap(parula); cb = colorbar;
ylabel(cb, 'NMSE_{in-band} (dB)');
caxis([-105, max(NMSE_ib_disp(:)) + 3]);
set(gca, 'XTick', 1:nM, 'XTickLabel', ...
    arrayfun(@(m) sprintf('%dX', m), M_vec, 'UniformOutput', false));
set(gca, 'YTick', 1:nK, 'YTickLabel', ...
    arrayfun(@(k) sprintf('K=%d', k), K_vec, 'UniformOutput', false));
xlabel('Oversampling  M = Fs / BW');
ylabel('DPD Polynomial Order  K');
title('In-Band Waveform Error (dB): only |f| < BW/2');

for ii = 1:nK
    for jj = 1:nM
        v = NMSE_inband(ii, jj);
        if v < -100, s = '0'; else, s = sprintf('%.1f', v); end
        if v > -60, clr = 'w'; else, clr = 'k'; end
        text(jj, ii, s, 'HorizontalAlignment', 'center', ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', clr);
        if K_vec(ii) > M_vec(jj)
            rectangle('Position', [jj-0.5 ii-0.5 1 1], ...
                      'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
        end
    end
end
saveas(gcf, 'DPD_sampling_NMSE_inband.png');
fprintf('Saved DPD_sampling_NMSE_inband.png\n');

%% 7. Figure 4 -- ACLR impact (pass basis through PA, compare ACLR)
ACLR1_grid = NaN(nK, nM);

for ii = 1:nK
    K = K_vec(ii);
    for jj = 1:nM
        M  = M_vec(jj);
        D  = D_vec(jj);
        Fs = M * BW;
        N_lo = N_ref / D;

        % PA output of basis at this Fs
        x_lo = x_ref(1:D:end);
        y_lo = x_lo .* abs(x_lo).^(K-1);
        y_pa = pa_poly(y_lo, pa_c);

        % PSD and ACLR at this Fs
        Pk = compute_psd(y_pa, N_lo);
        f_lo = ((-N_lo/2):(N_lo/2-1))' * (Fs / N_lo);
        ACLR1_grid(ii, jj) = compute_aclr1(Pk, f_lo, BW);
    end
end

figure('Position', [120 120 640 500]);
aclr_disp = ACLR1_grid;
aclr_disp(isnan(aclr_disp)) = 0;
imagesc(aclr_disp);
colormap(parula); cb = colorbar;
ylabel(cb, 'ACLR1 (dBc)');
set(gca, 'XTick', 1:nM, 'XTickLabel', ...
    arrayfun(@(m) sprintf('%dX', m), M_vec, 'UniformOutput', false));
set(gca, 'YTick', 1:nK, 'YTickLabel', ...
    arrayfun(@(k) sprintf('K=%d', k), K_vec, 'UniformOutput', false));
xlabel('Oversampling  M = Fs / BW');
ylabel('DPD Polynomial Order  K');
title('ACLR1 after PA — Basis at Fs');

for ii = 1:nK
    for jj = 1:nM
        v = ACLR1_grid(ii, jj);
        if isnan(v), s = 'N/A'; else, s = sprintf('%.1f', v); end
        if isnan(v) || v > median(ACLR1_grid(~isnan(ACLR1_grid(:))))
            clr = 'w';
        else
            clr = 'k';
        end
        text(jj, ii, s, 'HorizontalAlignment', 'center', ...
             'FontSize', 10, 'FontWeight', 'bold', 'Color', clr);
        if K_vec(ii) > M_vec(jj)
            rectangle('Position', [jj-0.5 ii-0.5 1 1], ...
                      'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
        end
    end
end

saveas(gcf, 'DPD_sampling_ACLR.png');
fprintf('Saved DPD_sampling_ACLR.png\n');

%% 8. Summary Tables
fprintf('\n=== Aliasing Condition: K > M ===\n');
fprintf('%8s', '');
for jj = 1:nM, fprintf(' %8s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    fprintf('K=%-5d', K_vec(ii));
    for jj = 1:nM
        if     K_vec(ii) < M_vec(jj), tag = 'OK';
        elseif K_vec(ii) == M_vec(jj), tag = 'Edge';
        else,                          tag = 'ALIASED';
        end
        fprintf(' %8s', tag);
    end
    fprintf('\n');
end

fprintf('\n=== In-Band NMSE (dB) — waveform error vs golden ===\n');
fprintf('%8s', '');
for jj = 1:nM, fprintf(' %10s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    fprintf('K=%-5d', K_vec(ii));
    for jj = 1:nM
        v = NMSE_inband(ii, jj);
        if v < -100, fprintf(' %10s', '<-100');
        else,        fprintf(' %10.2f', v);
        end
    end
    fprintf('\n');
end

fprintf('\n=== Total NMSE (dB) — full bandwidth ===\n');
fprintf('%8s', '');
for jj = 1:nM, fprintf(' %10s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    fprintf('K=%-5d', K_vec(ii));
    for jj = 1:nM
        v = NMSE_total(ii, jj);
        if v < -100, fprintf(' %10s', '<-100');
        else,        fprintf(' %10.2f', v);
        end
    end
    fprintf('\n');
end

fprintf('\nKey: K <= M -> no aliasing (NMSE ~ -inf)\n');
fprintf('     K >  M -> spectral content beyond Nyquist folds back\n');
fprintf('     Basis bandwidth = K * BW = K * %.0f MHz\n', BW/1e6);
fprintf('     Nyquist = Fs/2 = M * %.0f MHz\n', BW/2/1e6);

%% =========================================================================
%%  GAP ANALYSES (Sections 9-15)
%% =========================================================================

%% 9. Precise Aliasing Condition: K > M vs K >= 2M
fprintf('\n=== Precise Aliasing Condition ===\n');
fprintf('Basis bandwidth = K*BW.  Content at freq f folds to f-Fs.\n');
fprintf('In-band aliasing requires folded freq in [-BW/2, BW/2].\n');
fprintf('This happens when K*BW/2 > Fs - BW/2, i.e. K > 2M-1, i.e. K >= 2M.\n\n');

fprintf('%8s', '');
for jj = 1:nM, fprintf(' %10s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    K = K_vec(ii);
    fprintf('K=%-5d', K);
    for jj = 1:nM
        M = M_vec(jj);
        if K < M,       tag = 'OK';
        elseif K == M,   tag = 'Edge';
        elseif K < 2*M,  tag = 'OOB-only';
        else,            tag = 'IN-BAND';
        end
        fprintf(' %10s', tag);
    end
    fprintf('\n');
end
fprintf('\nOOB-only = aliasing lands in adjacent channel, not signal band\n');
fprintf('IN-BAND  = aliasing corrupts signal band (+-BW/2)\n');

%% 10. End-to-End DPD + PA Performance
fprintf('\n=== End-to-End DPD + PA ===\n');
n_iter   = 20;
G_target = 1.0;
z_cap    = 1.5;

DPD_NMSE  = NaN(nK, nM);
DPD_ACLR1 = NaN(nK, nM);
dpd_w     = cell(nK, nM);    % store coefficients

for jj = 1:nM
    D  = D_vec(jj);
    Fs = M_vec(jj) * BW;
    x_lo = x_ref(1:D:end);
    N_lo = length(x_lo);
    f_lo = ((-N_lo/2):(N_lo/2-1))' * (Fs / N_lo);

    for ii = 1:nK
        K  = K_vec(ii);
        A  = build_basis(x_lo, K);
        nc = size(A, 2);

        % Iterative multiplicative correction DPD
        w = zeros(nc, 1); w(1) = 1.0;
        for it = 1:n_iter
            z = A * w;
            y = pa_poly(z, pa_c);
            ratio = G_target * x_lo ./ (y + 1e-12 * (abs(y) < 1e-12));
            z_corr = z .* ratio;
            amp = abs(z_corr);
            over = amp > z_cap;
            z_corr(over) = z_corr(over) .* (z_cap ./ amp(over));
            w = (A' * A + lam * eye(nc)) \ (A' * z_corr);
        end
        dpd_w{ii, jj} = w;

        % Final evaluation
        z     = A * w;
        y_lin = pa_poly(z, pa_c);
        G_est = (y_lin' * x_lo) / (x_lo' * x_lo);
        err   = y_lin - G_est * x_lo;
        DPD_NMSE(ii, jj) = 10*log10(mean(abs(err).^2) / mean(abs(G_est*x_lo).^2));
        Pk = compute_psd(y_lin, N_lo);
        DPD_ACLR1(ii, jj) = compute_aclr1(Pk, f_lo, BW);
    end
end

% Print DPD NMSE table
fprintf('\n--- DPD Linearization NMSE (dB) ---\n');
fprintf('%8s', '');
for jj = 1:nM, fprintf(' %8s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    fprintf('K=%-5d', K_vec(ii));
    for jj = 1:nM, fprintf(' %8.2f', DPD_NMSE(ii,jj)); end
    fprintf('\n');
end

% DPD NMSE Heatmap
DPD_disp = max(DPD_NMSE, -60);
figure('Position', [130 130 640 500]);
imagesc(DPD_disp);
colormap(parula); cb = colorbar;
ylabel(cb, 'DPD NMSE (dB)');
caxis([min(DPD_disp(:))-2, max(DPD_disp(:))+2]);
set(gca, 'XTick', 1:nM, 'XTickLabel', ...
    arrayfun(@(m) sprintf('%dX', m), M_vec, 'UniformOutput', false));
set(gca, 'YTick', 1:nK, 'YTickLabel', ...
    arrayfun(@(k) sprintf('K=%d', k), K_vec, 'UniformOutput', false));
xlabel('Oversampling M'); ylabel('DPD Order K');
title('End-to-End DPD+PA NMSE (dB)');
for ii = 1:nK
    for jj = 1:nM
        v = DPD_NMSE(ii,jj);
        if v > median(DPD_NMSE(:)), clr='w'; else, clr='k'; end
        text(jj, ii, sprintf('%.1f', v), 'HorizontalAlignment','center', ...
             'FontSize', 10, 'FontWeight','bold', 'Color', clr);
        if K_vec(ii) >= 2*M_vec(jj)
            rectangle('Position',[jj-0.5 ii-0.5 1 1], ...
                      'EdgeColor','r','LineWidth',2);
        elseif K_vec(ii) > M_vec(jj)
            rectangle('Position',[jj-0.5 ii-0.5 1 1], ...
                      'EdgeColor',[1 0.5 0],'LineWidth',1.5,'LineStyle','--');
        end
    end
end
saveas(gcf, 'DPD_sampling_dpd_perf.png');
fprintf('Saved DPD_sampling_dpd_perf.png\n');

%% 11. DPD Coefficient Analysis
% For K=9 (max order), compare coefficients across M
fprintf('\n=== DPD Coefficients (K=9) ===\n');
idx_K9 = find(K_vec == 9);
odd_orders = 1:2:9;

figure('Position', [140 140 700 400]);
hold on;
colors_M = [0.85 0.20 0.20; 0.90 0.60 0.10; 0.20 0.60 0.80; 0.20 0.65 0.25];
for jj = 1:nM
    w9 = dpd_w{idx_K9, jj};
    w_mag = 20*log10(abs(w9) + 1e-15);
    plot(odd_orders, w_mag, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', colors_M(jj,:), 'MarkerFaceColor', colors_M(jj,:));
    fprintf('  M=%2d: |w| = [%s] (dB)\n', M_vec(jj), ...
            strjoin(arrayfun(@(v) sprintf('%.1f',v), w_mag, 'UniformOutput',false), ', '));
end
xlabel('Polynomial Order k'); ylabel('|w_k| (dB)');
title('DPD Coefficient Magnitudes (K=9 DPD)');
legend(arrayfun(@(m) sprintf('M=%d', m), M_vec, 'UniformOutput', false), ...
       'Location', 'southwest');
grid on;
set(gca, 'XTick', odd_orders);
saveas(gcf, 'DPD_sampling_coeff.png');
fprintf('Saved DPD_sampling_coeff.png\n');

%% 12. Observation Path Aliasing
fprintf('\n=== Observation Path Aliasing ===\n');
fprintf('PA output has BW up to 7*BW=140 MHz.  Feedback ADC at Fs aliases it.\n\n');

y_pa_ref = pa_poly(x_ref, pa_c);
Y_pa_ref = fftshift(fft(y_pa_ref)) / N_ref;

obs_NMSE_total  = NaN(1, nM);
obs_NMSE_inband = NaN(1, nM);

for jj = 1:nM
    D = D_vec(jj);
    y_pa_lo = y_pa_ref(1:D:end);
    y_pa_up = interpft(y_pa_lo, N_ref);
    Y_pa_up = fftshift(fft(y_pa_up)) / N_ref;

    err = y_pa_up - y_pa_ref;
    obs_NMSE_total(jj) = 10*log10(mean(abs(err).^2) / mean(abs(y_pa_ref).^2));

    err_ib = sum(abs(Y_pa_up(ib_mask) - Y_pa_ref(ib_mask)).^2);
    ref_ib = sum(abs(Y_pa_ref(ib_mask)).^2);
    obs_NMSE_inband(jj) = 10*log10(err_ib / ref_ib);

    fprintf('  Fs=%3dMHz (M=%2d): obs NMSE_total=%.2f dB, obs NMSE_inband=%.2f dB\n', ...
            M_vec(jj)*BW/1e6, M_vec(jj), obs_NMSE_total(jj), obs_NMSE_inband(jj));
end

%% 13. Basis Matrix Condition Number
fprintf('\n=== Basis Matrix Condition Number ===\n');
cond_grid = NaN(nK, nM);

for jj = 1:nM
    D = D_vec(jj);
    x_lo = x_ref(1:D:end);
    for ii = 1:nK
        K = K_vec(ii);
        A = build_basis(x_lo, K);
        cond_grid(ii, jj) = cond(A);
    end
end

fprintf('%8s', '');
for jj = 1:nM, fprintf(' %10s', sprintf('M=%d', M_vec(jj))); end
fprintf('\n');
for ii = 1:nK
    fprintf('K=%-5d', K_vec(ii));
    for jj = 1:nM, fprintf(' %10.1f', cond_grid(ii,jj)); end
    fprintf('\n');
end

figure('Position', [150 150 640 500]);
imagesc(log10(cond_grid));
colormap(hot); cb = colorbar;
ylabel(cb, 'log_{10}(cond)');
set(gca, 'XTick', 1:nM, 'XTickLabel', ...
    arrayfun(@(m) sprintf('%dX', m), M_vec, 'UniformOutput', false));
set(gca, 'YTick', 1:nK, 'YTickLabel', ...
    arrayfun(@(k) sprintf('K=%d', k), K_vec, 'UniformOutput', false));
xlabel('Oversampling M'); ylabel('DPD Order K');
title('Basis Matrix Condition Number');
for ii = 1:nK
    for jj = 1:nM
        v = cond_grid(ii,jj);
        if log10(v) > 2, clr='w'; else, clr='k'; end
        text(jj, ii, sprintf('%.0f', v), 'HorizontalAlignment','center', ...
             'FontSize', 9, 'FontWeight','bold', 'Color', clr);
    end
end
saveas(gcf, 'DPD_sampling_cond.png');
fprintf('Saved DPD_sampling_cond.png\n');

%% 14. K=9 @ 3X Deep Dive
fprintf('\n=== K=9 @ M=3 (Fs=60 MHz) Deep Dive ===\n');
K_dd = 9;  M_dd = 3;  Fs_dd = M_dd * BW;
D_dd = Fs_ref / Fs_dd;

x_dd  = x_ref(1:D_dd:end);
N_dd  = length(x_dd);
f_dd  = ((-N_dd/2):(N_dd/2-1))' * (Fs_dd / N_dd);

% Golden basis at Fs_ref
y_dd_ref = x_ref .* abs(x_ref).^(K_dd - 1);
Y_dd_ref = fftshift(fft(y_dd_ref)) / N_ref;
psd_dd_ref = abs(Y_dd_ref).^2;

% Aliased basis at 3X
y_dd_lo = x_dd .* abs(x_dd).^(K_dd - 1);
y_dd_up = interpft(y_dd_lo, N_ref);
Y_dd_up = fftshift(fft(y_dd_up)) / N_ref;
psd_dd_up = abs(Y_dd_up).^2;

% Error spectrum
psd_dd_err = abs(Y_dd_up - Y_dd_ref).^2;

% --- Figure: Spectral folding diagram ---
ref_pk = max(psd_dd_ref);

figure('Position', [40 40 1400 850]);

% 14a: Golden spectrum with band labels
subplot(2,2,1);
sdb = 10*log10(psd_dd_ref / ref_pk + 1e-30);
plot(f_ref_MHz, sdb, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2); hold on;

% Color bands by folding destination
bands = [0 10; 10 30; 30 50; 50 70; 70 90];
band_colors = [0.2 0.7 0.2;   % in-band (green)
               0.2 0.5 0.9;   % ACLR (blue)
               0.9 0.6 0.1;   % folds to ACLR (orange)
               0.85 0.15 0.15; % folds to IN-BAND (red!)
               0.6 0.3 0.7];  % folds to ACLR (purple)
band_labels = {'in-band (stays)', 'ACLR1 (stays)', ...
               '→ ACLR neg side', '→ IN-BAND!', '→ ACLR1'};
for b = 1:5
    % positive side
    idx = f_ref_MHz >= bands(b,1) & f_ref_MHz <= bands(b,2);
    area(f_ref_MHz(idx), max(sdb(idx), -80), 'BaseValue', -80, ...
         'FaceColor', band_colors(b,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    % negative side
    idx = f_ref_MHz >= -bands(b,2) & f_ref_MHz <= -bands(b,1);
    area(f_ref_MHz(idx), max(sdb(idx), -80), 'BaseValue', -80, ...
         'FaceColor', band_colors(b,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
xline(30, 'r--', 'Nyq', 'LineWidth', 1.5);
xline(-30, 'r--', 'LineWidth', 1.5);
xlim([-100 100]); ylim([-80 5]); grid on;
xlabel('MHz'); ylabel('PSD (dB)');
title('K=9 Golden Spectrum + Folding Bands');
legend(band_labels, 'Location', 'northeast', 'FontSize', 7);

% 14b: Result at 3X (aliased + reconstructed)
subplot(2,2,2);
sdb_up = 10*log10(psd_dd_up / ref_pk + 1e-30);
sdb_ref_clip = 10*log10(psd_dd_ref / ref_pk + 1e-30);
plot(f_ref_MHz, sdb_ref_clip, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); hold on;
plot(f_ref_MHz, sdb_up, 'b-', 'LineWidth', 1.2);
xline(30, 'r--', 'LineWidth', 1); xline(-30, 'r--', 'LineWidth', 1);
xline(10, 'k:', 'LineWidth', 0.8); xline(-10, 'k:', 'LineWidth', 0.8);
xlim([-40 40]); ylim([-60 5]); grid on;
xlabel('MHz'); ylabel('PSD (dB)');
title('After Decimation+Reconstruction (blue) vs Golden (gray)');
legend('Golden', 'At 3X (aliased)', 'Location', 'southwest');

% 14c: In-band zoom
subplot(2,2,3);
plot(f_ref_MHz, sdb_ref_clip, 'Color', [0.5 0.5 0.5], 'LineWidth', 1); hold on;
plot(f_ref_MHz, sdb_up, 'b-', 'LineWidth', 1.5);
sdb_err = 10*log10(psd_dd_err / ref_pk + 1e-30);
plot(f_ref_MHz, sdb_err, 'r-', 'LineWidth', 1);
xline(10, 'k--', 'LineWidth', 0.8); xline(-10, 'k--', 'LineWidth', 0.8);
xlim([-15 15]); ylim([-60 5]); grid on;
xlabel('MHz'); ylabel('PSD (dB)');
title('In-Band Zoom: Error from [50-70] MHz fold-back');
legend('Golden', 'Aliased', 'Error', 'Location', 'south');

% 14d: DPD comparison — K=9 at 3X vs 12X
subplot(2,2,4);
% DPD at 3X
w_3x = dpd_w{idx_K9, M_vec==3};
A_3x = build_basis(x_dd, 9);
z_3x = A_3x * w_3x;
y_3x = pa_poly(z_3x, pa_c);
Pk_3x = compute_psd(y_3x, N_dd);
sdb_3x = 10*log10(Pk_3x / max(Pk_3x) + 1e-30);

% DPD at 12X
x_12 = x_ref(1:D_vec(M_vec==12):end);
N_12 = length(x_12);
f_12 = ((-N_12/2):(N_12/2-1))' * (12*BW/N_12);
w_12 = dpd_w{idx_K9, M_vec==12};
A_12 = build_basis(x_12, 9);
z_12 = A_12 * w_12;
y_12 = pa_poly(z_12, pa_c);
Pk_12 = compute_psd(y_12, N_12);
sdb_12 = 10*log10(Pk_12 / max(Pk_12) + 1e-30);

% PA-only at 12X
y_pa_12 = pa_poly(x_12, pa_c);
Pk_pa12 = compute_psd(y_pa_12, N_12);
sdb_pa12 = 10*log10(Pk_pa12 / max(Pk_pa12) + 1e-30);

plot(f_12/1e6, sdb_pa12, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); hold on;
plot(f_12/1e6, sdb_12, 'Color', [0.2 0.65 0.3], 'LineWidth', 1.5);
plot(f_dd/1e6, sdb_3x, 'r-', 'LineWidth', 1.5);
xline(10, 'k:', 'LineWidth', 0.8); xline(-10, 'k:', 'LineWidth', 0.8);
xlim([-35 35]); ylim([-60 5]); grid on;
xlabel('MHz'); ylabel('PSD (dB, normalized)');
title(sprintf('K=9 DPD: 3X (NMSE=%.1f) vs 12X (NMSE=%.1f)', ...
              DPD_NMSE(idx_K9, M_vec==3), DPD_NMSE(idx_K9, M_vec==12)));
legend('PA only', 'DPD@12X', 'DPD@3X', 'Location', 'southwest');

sgtitle('K=9 @ M=3: Spectral Folding & DPD Impact');
saveas(gcf, 'DPD_sampling_K9_fold.png');
fprintf('Saved DPD_sampling_K9_fold.png\n');

% Print folding table
fprintf('\n--- K=9 @ 3X Spectral Folding ---\n');
fprintf('  Original band      → Folds to           Destination\n');
fprintf('  [  0, 10] MHz      → [ 0, 10] MHz       Signal band (stays)\n');
fprintf('  [ 10, 30] MHz      → [10, 30] MHz       ACLR1 (stays)\n');
fprintf('  [ 30, 50] MHz      → [-30,-10] MHz      ACLR1 neg (alias)\n');
fprintf('  [ 50, 70] MHz      → [-10, 10] MHz      IN-BAND (alias!)\n');
fprintf('  [ 70, 90] MHz      → [ 10, 30] MHz      ACLR1 (alias)\n');
fprintf('  Symmetric for negative frequencies.\n');

%% 15. Final Comprehensive Summary
fprintf('\n========================================\n');
fprintf('  COMPREHENSIVE SUMMARY\n');
fprintf('========================================\n');

fprintf('\n[1] Basis Aliasing: K > M causes total waveform error.\n');
fprintf('    But in-band aliasing requires K >= 2M.\n');

fprintf('\n[2] In-band aliasing only affects:\n');
fprintf('    K=7 @ M=3: NMSE_ib = %.1f dB\n', NMSE_inband(K_vec==7, M_vec==3));
fprintf('    K=9 @ M=3: NMSE_ib = %.1f dB\n', NMSE_inband(K_vec==9, M_vec==3));
fprintf('    All other (K,M) combinations: in-band NMSE < -100 dB\n');

fprintf('\n[3] End-to-end DPD performance (best K per M):\n');
for jj = 1:nM
    [best_nmse, best_idx] = min(DPD_NMSE(:, jj));
    fprintf('    M=%2d: best K=%d (NMSE=%.1f dB)\n', M_vec(jj), K_vec(best_idx), best_nmse);
end

fprintf('\n[4] Observation path aliasing (PA output BW=140MHz):\n');
for jj = 1:nM
    fprintf('    M=%2d: obs NMSE_total=%.1f dB, obs NMSE_inband=%.1f dB\n', ...
            M_vec(jj), obs_NMSE_total(jj), obs_NMSE_inband(jj));
end

fprintf('\n[5] Design rules:\n');
fprintf('    - M >= 6 : safe for all K <= 9 (no in-band aliasing)\n');
fprintf('    - M = 3  : safe for K <= 5 (OOB aliasing only, tolerable)\n');
fprintf('    -          K=7/9 @ M=3 has in-band aliasing, avoid.\n');

fprintf('\nAll figures saved.  Done.\n');

%% ===================== Local Functions =====================

function y = pa_poly(x, c)
    r2 = abs(x).^2;
    y  = x .* (c(1) + c(2)*r2 + c(3)*r2.^2 + c(4)*r2.^3);
end

function A = build_basis(x, K)
    orders = 1:2:K;
    A = zeros(length(x), length(orders));
    for i = 1:length(orders)
        A(:,i) = x .* abs(x).^(orders(i) - 1);
    end
end

function Pk = compute_psd(x, N)
    w  = hann(N, 'periodic');
    w  = w / sqrt(mean(w.^2));
    Xw = fftshift(fft(x .* w)) / N;
    Pk = abs(Xw).^2;
end

function a1 = compute_aclr1(Pk, f, BW)
    ib  = abs(f) <= BW/2;
    ac1 = (abs(f) > BW/2) & (abs(f) <= 3*BW/2);
    Pin = sum(Pk(ib));
    P1  = sum(Pk(ac1));
    if ~any(ac1) || P1 < 1e-25, a1 = NaN; else, a1 = 10*log10(P1/Pin); end
end
