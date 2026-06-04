% =========================================================
% MODULE 5c (corrected) — ERD onset vs EMG onset latency
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% CORRECTIONS vs previous version:
%   1. ERD onset detection: uses beta power TIME SERIES
%      (spectrogram sliding window) smoothed with moving
%      average, threshold = baseline mean - 1*baseline SD
%      sustained for at least 3 consecutive windows
%      → more robust than single-sample threshold crossing
%   2. EMG onset: identical to M03/M02c (mean+2SD baseline)
%   3. Delay = EMG_onset - ERD_onset (positive = ERD first)
%   4. Only trials where ERD_onset < EMG_onset kept
%      (physiologically: cortex desynchronizes before muscle)
%
% Expected: if tDCS speeds corticomotor transmission,
% delay POST < delay PRE in tDCS, larger reduction than sham
%
% Input:
%   - *_EEG_ica.set (sets_clean)
%   - *_EMG.set (sets)
%
% Output:
%   - Delay distribution figures per session
%   - Summary comparison across sessions
%   - all_latency_v2.mat
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

% Baseline window for ERD threshold (same as M02b)
S_BASE_ERD = round(1.5 * FS) + 1;
E_BASE_ERD = round(2.5 * FS);

% Search window for ERD onset: from stimulus to end of execution
S_SEARCH_ERD = round(3.0 * FS) + 1;
E_SEARCH_ERD = round(7.5 * FS);

% Baseline window for EMG threshold (same as M03)
S_BASE_EMG = round(1.5 * FS) + 1;
E_BASE_EMG = round(2.5 * FS);

% Search window for EMG onset
S_SEARCH_EMG = round(3.0 * FS) + 1;
E_SEARCH_EMG = round(8.0 * FS);

% Spectrogram parameters for beta time series
% Same as M02b: 0.5s window, 90% overlap
WIN_SPEC = round(0.5 * FS);    % 1200 samples
OVL_SPEC = round(0.45 * FS);   % 90% overlap
NFFT     = 2048;

BETA = [13, 30];

% ERD onset detection parameters
ERD_SD_THRESHOLD  = 1.0;   % threshold = baseline_mean - N*baseline_SD
ERD_MIN_SUSTAINED = 3;     % minimum consecutive windows below threshold

% Motor channels: C3=1, Cz=3, C4=5
MOTOR_CH = [1, 3, 5];

% EMG filter — identical to M03
NYQ = FS / 2;
[B_BP, A_BP] = butter(4, [20 500]/NYQ, 'bandpass');
[B_LP, A_LP] = butter(4, 8/NYQ, 'low');

% ── Session table ──────────────────────────────────────────
session_info;

all_latency = struct();

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
    % STEP 1: ERD onset per trial
    % ══════════════════════════════════════════════════════
    fprintf('  [1] Detecting ERD onset per trial...\n');

    [erd_onset_pre,  n_erd_pre]  = detect_erd_onset_block(pre_eeg,  eeg_folder, ...
        FS, TRIAL_SAMPLES, S_BASE_ERD, E_BASE_ERD, S_SEARCH_ERD, E_SEARCH_ERD, ...
        WIN_SPEC, OVL_SPEC, NFFT, BETA, MOTOR_CH, ERD_SD_THRESHOLD, ERD_MIN_SUSTAINED);

    [erd_onset_post, n_erd_post] = detect_erd_onset_block(post_eeg, eeg_folder, ...
        FS, TRIAL_SAMPLES, S_BASE_ERD, E_BASE_ERD, S_SEARCH_ERD, E_SEARCH_ERD, ...
        WIN_SPEC, OVL_SPEC, NFFT, BETA, MOTOR_CH, ERD_SD_THRESHOLD, ERD_MIN_SUSTAINED);

    fprintf('    ERD onset PRE:  %d/%d valid  mean=%.2fs  std=%.2fs\n', ...
        n_erd_pre,  length(erd_onset_pre),  mean(erd_onset_pre,'omitnan'), std(erd_onset_pre,'omitnan'));
    fprintf('    ERD onset POST: %d/%d valid  mean=%.2fs  std=%.2fs\n', ...
        n_erd_post, length(erd_onset_post), mean(erd_onset_post,'omitnan'), std(erd_onset_post,'omitnan'));

    % ══════════════════════════════════════════════════════
    % STEP 2: EMG onset per trial — identical to M03
    % ══════════════════════════════════════════════════════
    fprintf('  [2] Detecting EMG onset per trial...\n');

    [emg_onset_pre,  n_emg_pre]  = detect_emg_onset_block(pre_emg,  emg_folder, ...
        FS, TRIAL_SAMPLES, S_BASE_EMG, E_BASE_EMG, S_SEARCH_EMG, E_SEARCH_EMG, ...
        B_BP, A_BP, B_LP, A_LP);

    [emg_onset_post, n_emg_post] = detect_emg_onset_block(post_emg, emg_folder, ...
        FS, TRIAL_SAMPLES, S_BASE_EMG, E_BASE_EMG, S_SEARCH_EMG, E_SEARCH_EMG, ...
        B_BP, A_BP, B_LP, A_LP);

    fprintf('    EMG onset PRE:  %d/%d valid  mean=%.2fs\n', ...
        n_emg_pre,  length(emg_onset_pre),  mean(emg_onset_pre,'omitnan'));
    fprintf('    EMG onset POST: %d/%d valid  mean=%.2fs\n', ...
        n_emg_post, length(emg_onset_post), mean(emg_onset_post,'omitnan'));

    % ══════════════════════════════════════════════════════
    % STEP 3: Corticomotor delay = EMG onset - ERD onset
    % ══════════════════════════════════════════════════════
    fprintf('  [3] Computing corticomotor delay...\n');

    delay_pre  = emg_onset_pre  - erd_onset_pre;
    delay_post = emg_onset_post - erd_onset_post;

    % Keep only physiologically valid delays:
    % ERD must precede EMG (delay > 0) and be < 3s
    delay_pre(delay_pre  <= 0 | delay_pre  > 3) = NaN;
    delay_post(delay_post <= 0 | delay_post > 3) = NaN;

    n_delay_pre  = sum(~isnan(delay_pre));
    n_delay_post = sum(~isnan(delay_post));
    m_pre  = mean(delay_pre,  'omitnan');
    m_post = mean(delay_post, 'omitnan');
    s_pre  = std(delay_pre,   'omitnan');
    s_post = std(delay_post,  'omitnan');

    fprintf('    Delay PRE:  mean=%.3fs  std=%.3fs  n=%d/%d valid\n', ...
        m_pre,  s_pre,  n_delay_pre,  length(delay_pre));
    fprintf('    Delay POST: mean=%.3fs  std=%.3fs  n=%d/%d valid\n', ...
        m_post, s_post, n_delay_post, length(delay_post));
    fprintf('    Delta delay: %+.3fs  (%+.1f%%)\n', ...
        m_post - m_pre, 100*(m_post-m_pre)/m_pre);

    % Store
    all_latency.(label).erd_onset_pre   = erd_onset_pre;
    all_latency.(label).erd_onset_post  = erd_onset_post;
    all_latency.(label).emg_onset_pre   = emg_onset_pre;
    all_latency.(label).emg_onset_post  = emg_onset_post;
    all_latency.(label).delay_pre       = delay_pre;
    all_latency.(label).delay_post      = delay_post;
    all_latency.(label).mean_delay_pre  = m_pre;
    all_latency.(label).mean_delay_post = m_post;
    all_latency.(label).std_delay_pre   = s_pre;
    all_latency.(label).std_delay_post  = s_post;
    all_latency.(label).n_delay_pre     = n_delay_pre;
    all_latency.(label).n_delay_post    = n_delay_post;

    % ── Figures ───────────────────────────────────────────
    fig = plot_delay_summary(delay_pre, delay_post, m_pre, m_post, ...
        s_pre, s_post, n_delay_pre, n_delay_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_corticomotor_delay.png', label)));
    close(fig);

    fprintf('  [OK] Figures saved for %s\n', label);
end

% Summary across sessions
fig = plot_delay_comparison(all_latency);
saveas(fig, fullfile(ROOT_RESULTS, 'ALL_corticomotor_delay_summary.png'));
close(fig);

save(fullfile(ROOT_RESULTS, 'all_latency_v2.mat'), 'all_latency');
fprintf('\n=== Module 5c complete ===\n');


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


function [onsets, n_valid] = detect_erd_onset_block(file_list, folder, ...
        fs, trial_samples, s_base, e_base, s_search, e_search, ...
        win_spec, ovl_spec, nfft, beta, motor_ch, n_sd, min_sustained)
% Detects ERD onset per trial using beta power time series
%
% Method:
%   1. Compute spectrogram over full trial (same as M02b)
%   2. Extract beta power time series averaged across motor channels
%   3. Compute baseline: mean and SD of beta power in baseline window
%   4. Threshold = baseline_mean - n_sd * baseline_SD
%   5. Find first time point where beta stays below threshold
%      for at least min_sustained consecutive windows
%
% Returns onset in seconds within trial (NaN if not detected)

    onsets = [];

    for f = 1:length(file_list)
        EEG = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw = EEG.data;

        n_trials = floor(size(raw,2) / trial_samples);
        raw_ep   = reshape(raw(:,1:n_trials*trial_samples), ...
                           size(raw,1), trial_samples, n_trials);

        run_onsets = nan(n_trials, 1);

        for tr = 1:n_trials
            % Compute beta power time series averaged across motor channels
            beta_ts = [];
            t_ts    = [];

            for ci = 1:length(motor_ch)
                ch  = motor_ch(ci);
                sig = double(raw_ep(ch, :, tr));

                [~, f_ax, t_ax, P] = spectrogram(sig, hann(win_spec), ...
                    ovl_spec, nfft, fs);

                idx_beta = f_ax >= beta(1) & f_ax <= beta(2);
                bp = mean(P(idx_beta, :), 1);  % (1 x t_bins)

                if ci == 1
                    beta_ts = bp;
                    t_ts    = t_ax(:)';
                else
                    beta_ts = beta_ts + bp;
                end
            end
            beta_ts = beta_ts / length(motor_ch);

            % Baseline: time bins within baseline window
            idx_base = t_ts >= s_base/fs & t_ts <= e_base/fs;
            if sum(idx_base) < 2; continue; end

            base_mean = mean(beta_ts(idx_base));
            base_std  = std(beta_ts(idx_base));
            threshold = base_mean - n_sd * base_std;

            % Search window: time bins within search window
            idx_search = t_ts >= s_search/fs & t_ts <= e_search/fs;
            if sum(idx_search) < min_sustained; continue; end

            search_ts = beta_ts(idx_search);
            search_t  = t_ts(idx_search);

            % Find first sustained crossing below threshold
            below = search_ts < threshold;
            onset_idx = find_sustained_onset(below, min_sustained);

            if ~isnan(onset_idx)
                run_onsets(tr) = search_t(onset_idx);
            end
        end

        onsets = [onsets; run_onsets]; %#ok<AGROW>
    end

    n_valid = sum(~isnan(onsets));
end


function idx = find_sustained_onset(below_thresh, min_sustained)
% Returns index of first sustained crossing of min_sustained length
% Returns NaN if not found

    idx = NaN;
    n = length(below_thresh);

    for i = 1:(n - min_sustained + 1)
        if all(below_thresh(i:i+min_sustained-1))
            idx = i;
            return;
        end
    end
end


function [onsets, n_valid] = detect_emg_onset_block(file_list, folder, ...
        fs, trial_samples, s_base, e_base, s_search, e_search, ...
        b_bp, a_bp, b_lp, a_lp)
% EMG onset detection — identical to M03
% Returns onset in seconds within trial

    onsets = [];

    for f = 1:length(file_list)
        EEG = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw = EEG.data;

        % Bipolarize
        bipolar = zeros(4, size(raw,2));
        for m = 1:4
            bipolar(m,:) = raw(2*m-1,:) - raw(2*m,:);
        end

        % Bandpass + envelope
        for m = 1:4
            bipolar(m,:) = filtfilt(b_bp, a_bp, double(bipolar(m,:)));
            bipolar(m,:) = filtfilt(b_lp, a_lp, abs(bipolar(m,:)));
        end

        n_trials = floor(size(bipolar,2) / trial_samples);
        bp_ep    = reshape(bipolar(:,1:n_trials*trial_samples), ...
                           4, trial_samples, n_trials);

        run_onsets = nan(n_trials, 1);

        for tr = 1:n_trials
            env_mean  = mean(bp_ep(:,:,tr), 1);
            base_seg  = env_mean(s_base:e_base);
            threshold = mean(base_seg) + 2*std(base_seg);
            seg       = env_mean(s_search:e_search);
            cross     = find(seg > threshold);
            if ~isempty(cross)
                run_onsets(tr) = (s_search + cross(1) - 1) / fs;
            end
        end

        onsets = [onsets; run_onsets]; %#ok<AGROW>
    end

    n_valid = sum(~isnan(onsets));
end


% ── Plot functions ──────────────────────────────────────────

function fig = plot_delay_summary(delay_pre, delay_post, ...
        m_pre, m_post, s_pre, s_post, n_pre, n_post, label)

    fig = figure('Visible','off','Position',[100 100 1000 420]);

    subplot(1,2,1);
    valid_pre  = delay_pre(~isnan(delay_pre));
    valid_post = delay_post(~isnan(delay_post));

    if isempty(valid_pre) && isempty(valid_post)
        text(0.5,0.5,'No valid delays','HorizontalAlignment','center');
        axis off;
    else
        all_vals = [valid_pre; valid_post];
        edges = linspace(0, min(3, max(all_vals)*1.1), 25);
        histogram(valid_pre,  edges, 'FaceColor',[0.2 0.4 0.8], ...
            'FaceAlpha',0.6, 'DisplayName',sprintf('PRE (n=%d)',n_pre));
        hold on;
        histogram(valid_post, edges, 'FaceColor',[0.8 0.2 0.2], ...
            'FaceAlpha',0.6, 'DisplayName',sprintf('POST (n=%d)',n_post));
        if ~isnan(m_pre);  xline(m_pre,  'b--','LineWidth',2,'HandleVisibility','off'); end
        if ~isnan(m_post); xline(m_post, 'r--','LineWidth',2,'HandleVisibility','off'); end
        xlabel('Corticomotor delay (s)'); ylabel('Trial count');
        title('Delay distribution (ERD → EMG)');
        legend('Location','best'); grid on;
    end

    subplot(1,2,2);
    means = [m_pre, m_post];
    stds  = [s_pre, s_post];

    if any(isnan(means))
        text(0.5,0.5,'Insufficient data','HorizontalAlignment','center');
        axis off;
    else
        b = bar([1,2], means, 0.5);
        b.FaceColor = 'flat';
        b.CData = [0.2 0.4 0.8; 0.8 0.2 0.2];
        hold on;
        errorbar([1,2], means, stds, 'k', 'LineStyle','none', ...
            'LineWidth',1.5, 'CapSize',8);
        delta = m_post - m_pre;
        delta_pct = 100*delta/m_pre;
        text(1.5, max(means+stds)*1.15, ...
            sprintf('Δ=%+.3fs (%+.1f%%)', delta, delta_pct), ...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
        set(gca,'XTick',[1,2],'XTickLabel',{sprintf('PRE\n(n=%d)',n_pre), ...
            sprintf('POST\n(n=%d)',n_post)});
        ylabel('Mean corticomotor delay (s)');
        title('Mean ± SD'); grid on;
    end

    sgtitle(sprintf('Corticomotor delay (ERD→EMG) — %s', strrep(label,'_',' ')));
end


function fig = plot_delay_comparison(all_latency)

    fields = fieldnames(all_latency);
    n = length(fields);

    m_pre  = zeros(n,1); s_pre  = zeros(n,1); n_pre  = zeros(n,1);
    m_post = zeros(n,1); s_post = zeros(n,1); n_post = zeros(n,1);
    lbls   = cell(n,1);

    for i = 1:n
        m_pre(i)  = all_latency.(fields{i}).mean_delay_pre;
        m_post(i) = all_latency.(fields{i}).mean_delay_post;
        s_pre(i)  = all_latency.(fields{i}).std_delay_pre;
        s_post(i) = all_latency.(fields{i}).std_delay_post;
        n_pre(i)  = all_latency.(fields{i}).n_delay_pre;
        n_post(i) = all_latency.(fields{i}).n_delay_post;
        lbls{i}   = strrep(fields{i},'_',' ');
    end

    fig = figure('Visible','off','Position',[100 100 1000 450]);
    x = 1:n;

    bar(x-0.2, m_pre,  0.35, 'FaceColor',[0.2 0.4 0.8],'DisplayName','PRE');
    hold on;
    bar(x+0.2, m_post, 0.35, 'FaceColor',[0.8 0.2 0.2],'DisplayName','POST');
    errorbar(x-0.2, m_pre,  s_pre,  'k','LineStyle','none','LineWidth',1.2,'CapSize',6);
    errorbar(x+0.2, m_post, s_post, 'k','LineStyle','none','LineWidth',1.2,'CapSize',6);

    % Delta annotations
    for i = 1:n
        if ~isnan(m_pre(i)) && ~isnan(m_post(i))
            delta = m_post(i) - m_pre(i);
            y_top = max(m_pre(i)+s_pre(i), m_post(i)+s_post(i)) + 0.05;
            text(i, y_top, sprintf('%+.3fs', delta), ...
                'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
        end
    end

    set(gca,'XTick',x,'XTickLabel',lbls,'FontSize',9);
    ylabel('Mean corticomotor delay (s)');
    title('Corticomotor delay PRE vs POST — all sessions');
    legend('Location','best'); grid on;
end
