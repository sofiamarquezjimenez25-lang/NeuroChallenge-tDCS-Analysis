% =========================================================
% MODULE 5b (corrected) — ERD beta vs CCI correlation
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% CORRECTIONS vs previous version:
%   1. ERD per trial: pwelch with 1s window (not 0.5s)
%      → better freq resolution, more stable beta estimate
%   2. CCI per trial: exact same pipeline as M03
%      (same filter, same envelope, same RMS formula)
%      → consistent with reported CCI values
%
% For each trial computes:
%   ERD_beta(%) = (P_exec - P_base) / P_base * 100
%   CCI(%)      = (RMS_min / RMS_max) * 100
%
% Correlates ERD beta vs CCI across trials (Spearman).
% Question: do trials with deeper ERD show less cocontraction?
% Does this relationship change between tDCS and sham?
%
% Input:
%   - *_EEG_ica.set (sets_clean)
%   - *_EMG.set (sets)
%
% Output:
%   - Scatter plots ERD vs CCI PRE and POST
%   - Bar chart r values across sessions
%   - all_erd_cci_v2.mat
% =========================================================

clear; clc; close all;
eeglab nogui;

% ── Paths ──────────────────────────────────────────────────
ROOT_CLEAN   = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets_clean';
ROOT_SETS    = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets';
ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';

if ~isfolder(ROOT_RESULTS); mkdir(ROOT_RESULTS); end

% ── Constants ──────────────────────────────────────────────
FS            = 2400;
TRIAL_SAMPLES = FS * 10;

% Windows (samples within trial) — same as M02b and M03
S_BASE = round(1.5 * FS) + 1;   % baseline 1.5–2.5s
E_BASE = round(2.5 * FS);
S_EXEC = round(5.0 * FS) + 1;   % execution 5.0–7.5s
E_EXEC = round(7.5 * FS);

% pwelch parameters for ERD per trial
% 1s window over 2.5s execution → ~4 segments, ~1Hz freq resolution
WIN_PWELCH = round(1.0 * FS);   % 1s = 2400 samples
OVL_PWELCH = round(0.5 * FS);   % 50% overlap
NFFT       = 4096;               % fine freq resolution

BETA  = [13, 30];
ALPHA = [8,  13];

% Motor channels for ERD: C3=1, Cz=3, C4=5
MOTOR_CH = [1, 3, 5];

% EMG filter parameters — identical to M03
NYQ = FS / 2;
[B_BP, A_BP] = butter(4, [20 500]/NYQ, 'bandpass');
[B_LP, A_LP] = butter(4, 8/NYQ, 'low');

% ── Session table ──────────────────────────────────────────
session_info;

all_erd_cci = struct();

% ── Main loop ──────────────────────────────────────────────
for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);

    eeg_folder = fullfile(ROOT_CLEAN, participant, session);
    emg_folder = fullfile(ROOT_SETS,  participant, session);

    fprintf('\n=== %s ===\n', label);

    % ── Classify files ────────────────────────────────────
    ica_files     = dir(fullfile(eeg_folder, '*_EEG_ica.set'));
    emg_files_all = dir(fullfile(emg_folder, '*_EMG.set'));
    [~,idx] = sort({ica_files.name});     ica_files     = ica_files(idx);
    [~,idx] = sort({emg_files_all.name}); emg_files_all = emg_files_all(idx);

    pre_eeg = {}; post_eeg = {};
    pre_emg = {}; post_emg = {};

    for f = 1:length(ica_files)
        rn = get_run_number(ica_files(f).name);
        switch rn
            case {3,4,5,6},  pre_eeg{end+1}  = ica_files(f).name;
            case {7,8,9,10}, post_eeg{end+1} = ica_files(f).name;
        end
    end
    for f = 1:length(emg_files_all)
        rn = get_run_number(emg_files_all(f).name);
        switch rn
            case {3,4,5,6},  pre_emg{end+1}  = emg_files_all(f).name;
            case {7,8,9,10}, post_emg{end+1} = emg_files_all(f).name;
        end
    end

    if isempty(pre_eeg) || isempty(post_eeg) || isempty(pre_emg) || isempty(post_emg)
        fprintf('  [WARNING] Missing files\n'); continue;
    end

    % ══════════════════════════════════════════════════════
    % STEP 1: ERD beta per trial
    % ══════════════════════════════════════════════════════
    fprintf('  [1] ERD beta per trial...\n');

    erd_pre  = compute_erd_per_trial(pre_eeg,  eeg_folder, FS, TRIAL_SAMPLES, ...
        S_BASE, E_BASE, S_EXEC, E_EXEC, WIN_PWELCH, OVL_PWELCH, NFFT, BETA, MOTOR_CH);
    erd_post = compute_erd_per_trial(post_eeg, eeg_folder, FS, TRIAL_SAMPLES, ...
        S_BASE, E_BASE, S_EXEC, E_EXEC, WIN_PWELCH, OVL_PWELCH, NFFT, BETA, MOTOR_CH);

    fprintf('    PRE:  n=%d  mean=%.1f%%  std=%.1f%%\n', ...
        sum(~isnan(erd_pre)),  mean(erd_pre,'omitnan'),  std(erd_pre,'omitnan'));
    fprintf('    POST: n=%d  mean=%.1f%%  std=%.1f%%\n', ...
        sum(~isnan(erd_post)), mean(erd_post,'omitnan'), std(erd_post,'omitnan'));

    % ══════════════════════════════════════════════════════
    % STEP 2: CCI per trial — identical pipeline to M03
    % ══════════════════════════════════════════════════════
    fprintf('  [2] CCI per trial (M03-identical pipeline)...\n');

    cci_pre  = compute_cci_per_trial(pre_emg,  emg_folder, FS, TRIAL_SAMPLES, ...
        S_EXEC, E_EXEC, B_BP, A_BP, B_LP, A_LP);
    cci_post = compute_cci_per_trial(post_emg, emg_folder, FS, TRIAL_SAMPLES, ...
        S_EXEC, E_EXEC, B_BP, A_BP, B_LP, A_LP);

    fprintf('    PRE:  n=%d  mean=%.1f%%  std=%.1f%%\n', ...
        sum(~isnan(cci_pre)),  mean(cci_pre,'omitnan'),  std(cci_pre,'omitnan'));
    fprintf('    POST: n=%d  mean=%.1f%%  std=%.1f%%\n', ...
        sum(~isnan(cci_post)), mean(cci_post,'omitnan'), std(cci_post,'omitnan'));

    % Sanity check: block mean should match M03 values
    fprintf('    [CHECK] CCI block means — PRE=%.1f%% POST=%.1f%% (should match M03)\n', ...
        mean(cci_pre,'omitnan'), mean(cci_post,'omitnan'));

    % ══════════════════════════════════════════════════════
    % STEP 3: Spearman correlation ERD vs CCI
    % ══════════════════════════════════════════════════════
    fprintf('  [3] Correlating ERD beta vs CCI...\n');

    [r_pre,  p_pre,  n_pre]  = spearman_valid(erd_pre,  cci_pre);
    [r_post, p_post, n_post] = spearman_valid(erd_post, cci_post);

    fprintf('    PRE:  r=%+.3f  p=%.3f  n=%d  %s\n', ...
        r_pre,  p_pre,  n_pre,  sig_label(p_pre));
    fprintf('    POST: r=%+.3f  p=%.3f  n=%d  %s\n', ...
        r_post, p_post, n_post, sig_label(p_post));
    fprintf('    Delta r = %+.3f\n', r_post - r_pre);

    % Store
    all_erd_cci.(label).erd_pre      = erd_pre;
    all_erd_cci.(label).erd_post     = erd_post;
    all_erd_cci.(label).cci_pre      = cci_pre;
    all_erd_cci.(label).cci_post     = cci_post;
    all_erd_cci.(label).r_pre        = r_pre;
    all_erd_cci.(label).r_post       = r_post;
    all_erd_cci.(label).p_pre        = p_pre;
    all_erd_cci.(label).p_post       = p_post;
    all_erd_cci.(label).n_pre        = n_pre;
    all_erd_cci.(label).n_post       = n_post;

    % ══════════════════════════════════════════════════════
    % STEP 4: Figures
    % ══════════════════════════════════════════════════════

    % Figure 1: Scatter PRE and POST with regression line
    fig = plot_scatter(erd_pre, cci_pre, erd_post, cci_post, ...
        r_pre, p_pre, n_pre, r_post, p_post, n_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_ERD_CCI_scatter.png', label)));
    close(fig);

    % Figure 2: ERD and CCI distributions PRE vs POST
    fig = plot_distributions(erd_pre, erd_post, cci_pre, cci_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_ERD_CCI_distributions.png', label)));
    close(fig);

    fprintf('  [OK] Figures saved for %s\n', label);
end

% Summary bar chart of r values across sessions
fig = plot_r_summary(all_erd_cci);
saveas(fig, fullfile(ROOT_RESULTS, 'ALL_ERD_CCI_r_summary.png'));
close(fig);

save(fullfile(ROOT_RESULTS, 'all_erd_cci_v2.mat'), 'all_erd_cci');
fprintf('\n=== Module 5b complete ===\n');


% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function run_num = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok); run_num = str2double(tok{1}{1});
    else
        tok = regexp(fname, 'run(\d+)', 'tokens', 'ignorecase');
        if ~isempty(tok); run_num = str2double(tok{1}{1});
        else; run_num = 0; end
    end
end


function erd = compute_erd_per_trial(file_list, folder, fs, trial_samples, ...
        s_base, e_base, s_exec, e_exec, win_pw, ovl_pw, nfft, beta, motor_ch)
% Computes beta ERD (%) per trial averaged across motor channels
% ERD(%) = (P_exec_beta - P_base_beta) / P_base_beta * 100
% Uses pwelch with 1s window for stable spectral estimate

    erd = [];

    for f = 1:length(file_list)
        EEG = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw = EEG.data;  % (8 x samples)

        n_trials = floor(size(raw,2) / trial_samples);
        raw_ep   = reshape(raw(:,1:n_trials*trial_samples), ...
                           size(raw,1), trial_samples, n_trials);

        erd_run = nan(n_trials, 1);

        for tr = 1:n_trials
            ch_erd = nan(length(motor_ch), 1);

            for ci = 1:length(motor_ch)
                ch = motor_ch(ci);

                % Baseline segment
                base_seg = double(raw_ep(ch, s_base:e_base, tr));
                [pxx_b, f_b] = pwelch(base_seg, hann(min(win_pw, length(base_seg))), ...
                    ovl_pw, nfft, fs);
                idx_b = f_b >= beta(1) & f_b <= beta(2);
                p_base = mean(pxx_b(idx_b));

                if p_base <= 0; continue; end

                % Execution segment
                exec_seg = double(raw_ep(ch, s_exec:e_exec, tr));
                [pxx_e, f_e] = pwelch(exec_seg, hann(min(win_pw, length(exec_seg))), ...
                    ovl_pw, nfft, fs);
                idx_e = f_e >= beta(1) & f_e <= beta(2);
                p_exec = mean(pxx_e(idx_e));

                ch_erd(ci) = (p_exec - p_base) / p_base * 100;
            end

            % Average across motor channels
            erd_run(tr) = mean(ch_erd, 'omitnan');
        end

        erd = [erd; erd_run]; %#ok<AGROW>
    end
end


function cci = compute_cci_per_trial(file_list, folder, fs, trial_samples, ...
        s_exec, e_exec, b_bp, a_bp, b_lp, a_lp)
% Computes CCI per trial using IDENTICAL pipeline to M03
% CCI = (RMS_min / RMS_max) * 100 across M1-M4

    cci = [];

    for f = 1:length(file_list)
        EEG = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw = EEG.data;

        % Bipolarize — identical to M03
        bipolar = zeros(4, size(raw,2));
        for m = 1:4
            bipolar(m,:) = raw(2*m-1,:) - raw(2*m,:);
        end

        % Bandpass 20-500 Hz — identical to M03
        for m = 1:4
            bipolar(m,:) = filtfilt(b_bp, a_bp, double(bipolar(m,:)));
        end

        % Envelope: rectify + LP 8Hz — identical to M03
        envelope = zeros(4, size(bipolar,2));
        for m = 1:4
            envelope(m,:) = filtfilt(b_lp, a_lp, abs(bipolar(m,:)));
        end

        % Epoch
        n_trials = floor(size(envelope,2) / trial_samples);
        env_ep   = reshape(envelope(:,1:n_trials*trial_samples), ...
                           4, trial_samples, n_trials);

        cci_run = nan(n_trials, 1);

        for tr = 1:n_trials
            rms_m = zeros(4,1);
            for m = 1:4
                seg    = env_ep(m, s_exec:e_exec, tr);
                rms_m(m) = sqrt(mean(seg.^2));
            end

            act = rms_m;  % (4x1) RMS per muscle
            pair1 = 2 * min(act(1), act(3)) / (act(1) + act(3) + 1e-10) * 100;
            pair2 = 2 * min(act(2), act(4)) / (act(2) + act(4) + 1e-10) * 100;
            cci_run(tr) = (pair1 + pair2) / 2;
        end

        cci = [cci; cci_run]; %#ok<AGROW>
    end
end


function [r, p, n] = spearman_valid(x, y)
% Spearman correlation ignoring NaN pairs
    valid = ~isnan(x) & ~isnan(y);
    n = sum(valid);
    if n < 5; r = NaN; p = NaN; return; end
    [r, p] = corr(x(valid), y(valid), 'Type', 'Spearman');
end


function lbl = sig_label(p)
    if isnan(p);    lbl = 'n/a';
    elseif p < 0.001; lbl = '***';
    elseif p < 0.01;  lbl = '**';
    elseif p < 0.05;  lbl = '*';
    else;             lbl = 'ns';
    end
end


% ── Plot functions ──────────────────────────────────────────

function fig = plot_scatter(erd_pre, cci_pre, erd_post, cci_post, ...
        r_pre, p_pre, n_pre, r_post, p_post, n_post, label)

    fig = figure('Visible','off','Position',[100 100 1000 450]);

    data_erd = {erd_pre,  erd_post};
    data_cci = {cci_pre,  cci_post};
    rs       = [r_pre,    r_post];
    ps       = [p_pre,    p_post];
    ns       = [n_pre,    n_post];
    titles   = {'PRE','POST'};
    cols     = {[0.2 0.4 0.8], [0.8 0.2 0.2]};

    for b = 1:2
        subplot(1,2,b);
        x = data_erd{b};
        y = data_cci{b};
        valid = ~isnan(x) & ~isnan(y);

        scatter(x(valid), y(valid), 20, cols{b}, 'filled', ...
            'MarkerFaceAlpha', 0.4); hold on;

        % Regression line
        if sum(valid) > 2
            cf = polyfit(x(valid), y(valid), 1);
            xr = linspace(min(x(valid)), max(x(valid)), 100);
            plot(xr, polyval(cf,xr), '-', 'Color', cols{b}, 'LineWidth', 2);
        end

        % Reference lines
        xline(0, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');

        xlabel('Beta ERD (%)');
        ylabel('CCI (%)');
        title(sprintf('%s — r=%+.3f, p=%.3f %s (n=%d)', ...
            titles{b}, rs(b), ps(b), sig_label(ps(b)), ns(b)));
        grid on;
    end

    sgtitle(sprintf('ERD beta vs CCI — Spearman trial-by-trial — %s', ...
        strrep(label,'_',' ')));
end


function fig = plot_distributions(erd_pre, erd_post, cci_pre, cci_post, label)
% Distribution plots: ERD and CCI PRE vs POST

    fig = figure('Visible','off','Position',[100 100 1000 450]);

    % ERD distribution
    subplot(1,2,1);
    valid_pre  = erd_pre(~isnan(erd_pre));
    valid_post = erd_post(~isnan(erd_post));
    edges_erd  = linspace(min([valid_pre;valid_post]), ...
                          max([valid_pre;valid_post]), 30);
    histogram(valid_pre,  edges_erd, 'FaceColor',[0.2 0.4 0.8], ...
        'FaceAlpha',0.6, 'DisplayName',sprintf('PRE (n=%d)',length(valid_pre)));
    hold on;
    histogram(valid_post, edges_erd, 'FaceColor',[0.8 0.2 0.2], ...
        'FaceAlpha',0.6, 'DisplayName',sprintf('POST (n=%d)',length(valid_post)));
    xline(mean(valid_pre),  'b--', 'LineWidth',2, 'HandleVisibility','off');
    xline(mean(valid_post), 'r--', 'LineWidth',2, 'HandleVisibility','off');
    xlabel('Beta ERD (%)'); ylabel('Trial count');
    title('Beta ERD distribution');
    legend('Location','best'); grid on;

    % CCI distribution
    subplot(1,2,2);
    valid_cci_pre  = cci_pre(~isnan(cci_pre));
    valid_cci_post = cci_post(~isnan(cci_post));
    edges_cci = linspace(0, 100, 30);
    histogram(valid_cci_pre,  edges_cci, 'FaceColor',[0.2 0.4 0.8], ...
        'FaceAlpha',0.6, 'DisplayName',sprintf('PRE (n=%d)',length(valid_cci_pre)));
    hold on;
    histogram(valid_cci_post, edges_cci, 'FaceColor',[0.8 0.2 0.2], ...
        'FaceAlpha',0.6, 'DisplayName',sprintf('POST (n=%d)',length(valid_cci_post)));
    xline(mean(valid_cci_pre),  'b--', 'LineWidth',2, 'HandleVisibility','off');
    xline(mean(valid_cci_post), 'r--', 'LineWidth',2, 'HandleVisibility','off');
    xlabel('CCI (%)'); ylabel('Trial count');
    title('CCI distribution');
    legend('Location','best'); grid on;

    sgtitle(sprintf('ERD and CCI distributions — %s', strrep(label,'_',' ')));
end


function fig = plot_r_summary(all_erd_cci)
% Bar chart: r PRE vs POST across all sessions

    fields = fieldnames(all_erd_cci);
    n = length(fields);

    r_pre  = zeros(n,1);
    r_post = zeros(n,1);
    p_pre  = zeros(n,1);
    p_post = zeros(n,1);
    lbls   = cell(n,1);

    for i = 1:n
        r_pre(i)  = all_erd_cci.(fields{i}).r_pre;
        r_post(i) = all_erd_cci.(fields{i}).r_post;
        p_pre(i)  = all_erd_cci.(fields{i}).p_pre;
        p_post(i) = all_erd_cci.(fields{i}).p_post;
        lbls{i}   = strrep(fields{i},'_',' ');
    end

    fig = figure('Visible','off','Position',[100 100 900 450]);
    x = 1:n;
    bar(x-0.2, r_pre,  0.35, 'FaceColor',[0.2 0.4 0.8],'DisplayName','PRE');
    hold on;
    bar(x+0.2, r_post, 0.35, 'FaceColor',[0.8 0.2 0.2],'DisplayName','POST');
    yline(0,'k-','LineWidth',1.2);

    % Significance markers
    for i = 1:n
        y_pre  = r_pre(i)  + sign(r_pre(i)) * 0.02;
        y_post = r_post(i) + sign(r_post(i)) * 0.02;
        text(i-0.2, y_pre,  sig_label(p_pre(i)),  'HorizontalAlignment','center','FontSize',10);
        text(i+0.2, y_post, sig_label(p_post(i)), 'HorizontalAlignment','center','FontSize',10);
    end

    set(gca,'XTick',x,'XTickLabel',lbls,'FontSize',9);
    ylabel('Spearman r (ERD beta vs CCI)');
    title('ERD–CCI correlation PRE vs POST — all sessions');
    legend('Location','best'); grid on;
    ylim([min([r_pre;r_post])*1.3 - 0.05, max([r_pre;r_post])*1.3 + 0.05]);
end
