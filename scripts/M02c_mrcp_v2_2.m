% =========================================================
% MODULE 2c — MRCP (Movement-Related Cortical Potential)
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Computes MRCP aligned to EMG onset per trial.
% EMG onset detected internally (no external file needed).
%
% Filter: LP 4 Hz applied per run before concatenation
% Baseline: -1.5 to -1.0 s relative to EMG onset
% Epoch: -1.5 to +1.0 s
% Rejection: |amplitude| > 150 µV
%
% Metrics:
%   RP amplitude: mean Cz from -0.5 to 0 s
%   MP amplitude: min  Cz from  0   to +0.1 s
%
% NOTE: BP slow ramp not expected in cued tasks.
%       MP is the primary metric for this paradigm.
% =========================================================

clear; clc; close all;
eeglab nogui;

ROOT_CLEAN   = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets_clean';
ROOT_SETS    = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets';
ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';
if ~isfolder(ROOT_RESULTS); mkdir(ROOT_RESULTS); end

FS            = 2400;
TRIAL_SAMPLES = FS * 10;

% MRCP epoch parameters
MRCP_PRE_S  = 1.5;    % s before EMG onset
MRCP_POST_S = 1.0;    % s after  EMG onset
BASE_WIN    = [-1.5, -1.0];  % baseline window (pre-movement rest)
LP_CUTOFF   = 4;             % Hz

% Metric windows
RP_WIN = [-0.5,  0.0];
MP_WIN = [ 0.0,  0.1];

% EMG onset detection parameters (same as M03)
WIN_BASE_EMG = [1.5, 2.5];  % baseline window within trial

CH_NAMES = {'C3','C1','Cz','C2','C4','CP1','CP2','Pz'};
CZ_IDX   = 3;

session_info;
all_mrcp = struct();

for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);

    eeg_folder = fullfile(ROOT_CLEAN, participant, session);
    emg_folder = fullfile(ROOT_SETS,  participant, session);

    fprintf('\n=== %s ===\n', label);

    % Classify EEG files
    ica_files = dir(fullfile(eeg_folder, '*_EEG_ica.set'));
    if isempty(ica_files)
        fprintf('  [SKIP] No ICA files found\n'); continue;
    end
    [~, idx] = sort({ica_files.name});
    ica_files = ica_files(idx);

    % Classify EMG files
    emg_files = dir(fullfile(emg_folder, '*_EMG.set'));
    [~, idx]  = sort({emg_files.name});
    emg_files = emg_files(idx);

    pre_eeg = {}; post_eeg = {};
    pre_emg = {}; post_emg = {};

    for f = 1:length(ica_files)
        rn = get_run_number(ica_files(f).name);
        switch rn
            case {3,4,5,6},  pre_eeg{end+1}  = ica_files(f).name;
            case {7,8,9,10}, post_eeg{end+1} = ica_files(f).name;
        end
    end
    for f = 1:length(emg_files)
        rn = get_run_number(emg_files(f).name);
        switch rn
            case {3,4,5,6},  pre_emg{end+1}  = emg_files(f).name;
            case {7,8,9,10}, post_emg{end+1} = emg_files(f).name;
        end
    end

    if isempty(pre_eeg) || isempty(post_eeg) || isempty(pre_emg) || isempty(post_emg)
        fprintf('  [WARNING] Missing PRE or POST files\n'); continue;
    end

    % ── Detect EMG onsets ──────────────────────────────────
    onsets_pre  = detect_emg_onsets(pre_emg,  emg_folder, FS, TRIAL_SAMPLES, WIN_BASE_EMG);
    onsets_post = detect_emg_onsets(post_emg, emg_folder, FS, TRIAL_SAMPLES, WIN_BASE_EMG);

    fprintf('  EMG onsets: PRE=%d/%d valid  POST=%d/%d valid\n', ...
        sum(~isnan(onsets_pre)),  length(onsets_pre), ...
        sum(~isnan(onsets_post)), length(onsets_post));

    % ── Compute MRCP ──────────────────────────────────────
    [mrcp_pre,  n_pre,  t_mrcp] = compute_mrcp_block(pre_eeg,  eeg_folder, onsets_pre,  FS, TRIAL_SAMPLES, MRCP_PRE_S, MRCP_POST_S, BASE_WIN, LP_CUTOFF);
    [mrcp_post, n_post, ~     ] = compute_mrcp_block(post_eeg, eeg_folder, onsets_post, FS, TRIAL_SAMPLES, MRCP_PRE_S, MRCP_POST_S, BASE_WIN, LP_CUTOFF);

    fprintf('  Valid epochs: PRE=%d  POST=%d\n', n_pre, n_post);

    % ── Extract metrics at Cz ──────────────────────────────
    rp_idx = t_mrcp >= RP_WIN(1) & t_mrcp <= RP_WIN(2);
    mp_idx = t_mrcp >= MP_WIN(1) & t_mrcp <= MP_WIN(2);

    rp_pre  = mean(mrcp_pre(CZ_IDX,  rp_idx));
    rp_post = mean(mrcp_post(CZ_IDX, rp_idx));
    mp_pre  = min(mrcp_pre(CZ_IDX,   mp_idx));
    mp_post = min(mrcp_post(CZ_IDX,  mp_idx));

    fprintf('  Cz — RP: PRE=%.3f µV  POST=%.3f µV  Δ=%+.3f\n', rp_pre, rp_post, rp_post-rp_pre);
    fprintf('  Cz — MP: PRE=%.3f µV  POST=%.3f µV  Δ=%+.3f\n', mp_pre, mp_post, mp_post-mp_pre);

    % ── Store ──────────────────────────────────────────────
    all_mrcp.(label).mrcp_pre    = mrcp_pre;
    all_mrcp.(label).mrcp_post   = mrcp_post;
    all_mrcp.(label).t_mrcp      = t_mrcp;
    all_mrcp.(label).n_valid_pre  = n_pre;
    all_mrcp.(label).n_valid_post = n_post;
    all_mrcp.(label).metrics_pre.rp_cz  = rp_pre;
    all_mrcp.(label).metrics_pre.mp_cz  = mp_pre;
    all_mrcp.(label).metrics_post.rp_cz = rp_post;
    all_mrcp.(label).metrics_post.mp_cz = mp_post;
    all_mrcp.(label).onsets_pre   = onsets_pre;
    all_mrcp.(label).onsets_post  = onsets_post;

    % ── Figures ────────────────────────────────────────────
    fig = plot_mrcp(mrcp_pre, mrcp_post, t_mrcp, label, CH_NAMES, RP_WIN, MP_WIN, rp_pre, rp_post, mp_pre, mp_post);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_MRCP.png', label)));
    close(fig);

    fig = plot_mrcp_grid(mrcp_pre, mrcp_post, t_mrcp, label, CH_NAMES);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_MRCP_channels.png', label)));
    close(fig);

    fprintf('  [OK] Figures saved for %s\n', label);
end

save(fullfile(ROOT_RESULTS, 'all_mrcp.mat'), 'all_mrcp');
fprintf('\n=== Module 2c complete ===\n');


% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function rn = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok); rn = str2double(tok{1}{1});
    else; rn = 0; end
end


function onsets = detect_emg_onsets(file_list, folder, fs, trial_samples, win_base)
% Detects EMG onset per trial — same method as M03
% Returns onset in seconds within trial (NaN if not detected)

    nyq = fs/2;
    [b_bp, a_bp] = butter(4, [20 500]/nyq, 'bandpass');
    [b_lp, a_lp] = butter(4, 8/nyq, 'low');

    s_base   = round(win_base(1)*fs) + 1;
    e_base   = round(win_base(2)*fs);
    s_search = round(3.0*fs) + 1;
    e_search = round(8.0*fs);

    onsets = [];
    for f = 1:length(file_list)
        EEG = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw = double(EEG.data);

        bipolar = zeros(4, size(raw,2));
        for m = 1:4
            bipolar(m,:) = raw(2*m-1,:) - raw(2*m,:);
            bipolar(m,:) = filtfilt(b_bp, a_bp, bipolar(m,:));
            bipolar(m,:) = filtfilt(b_lp, a_lp, abs(bipolar(m,:)));
        end

        n_trials = floor(size(bipolar,2) / trial_samples);
        bp_ep    = reshape(bipolar(:,1:n_trials*trial_samples), 4, trial_samples, n_trials);

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
end


function [mrcp_avg, n_valid, t_mrcp] = compute_mrcp_block(eeg_files, folder, onsets, fs, trial_samples, pre_s, post_s, base_win, lp_cutoff)
% Computes MRCP averaged across valid trials
% EEG is filtered LP per run before concatenation

    epoch_len = round((pre_s + post_s) * fs) + 1;
    t_mrcp    = linspace(-pre_s, post_s, epoch_len);
    idx_base  = t_mrcp >= base_win(1) & t_mrcp <= base_win(2);

    nyq = fs/2;
    [b_lp, a_lp] = butter(2, lp_cutoff/nyq, 'low');

    % Load and filter per run to avoid edge artifacts at run boundaries
    data_all = [];
    for f = 1:length(eeg_files)
        EEG    = pop_loadset('filename', eeg_files{f}, 'filepath', folder);
        data_f = double(EEG.data);
        for ch = 1:size(data_f,1)
            data_f(ch,:) = filtfilt(b_lp, a_lp, data_f(ch,:));
        end
        data_all = [data_all, data_f]; %#ok<AGROW>
    end

    n_ch      = size(data_all, 1);
    mrcp_sum  = zeros(n_ch, epoch_len);
    n_valid   = 0;

    for tr = 1:length(onsets)
        if isnan(onsets(tr)); continue; end

        % Onset sample in concatenated signal
        trial_start     = (tr-1) * trial_samples;
        onset_in_concat = trial_start + round(onsets(tr) * fs);

        ep_start = onset_in_concat - round(pre_s  * fs);
        ep_end   = onset_in_concat + round(post_s * fs);

        if ep_start < 1 || ep_end > size(data_all,2); continue; end

        epoch = data_all(:, ep_start:ep_end);
        if size(epoch,2) ~= epoch_len; continue; end

        % Baseline correction
        baseline = mean(epoch(:, idx_base), 2);
        epoch    = epoch - baseline;

        % Artifact rejection
        if max(abs(epoch(:))) > 150; continue; end

        mrcp_sum = mrcp_sum + epoch;
        n_valid  = n_valid + 1;
    end

    if n_valid > 0
        mrcp_avg = mrcp_sum / n_valid;
    else
        mrcp_avg = zeros(n_ch, epoch_len);
        fprintf('    [WARNING] No valid epochs\n');
    end
end


% ── Plot functions ──────────────────────────────────────────

function fig = plot_mrcp(mrcp_pre, mrcp_post, t, label, ch_names, rp_win, mp_win, rp_pre, rp_post, mp_pre, mp_post)
% Main MRCP figure: Cz PRE vs POST with metric annotations

    cz = 3;
    fig = figure('Visible','off','Position',[100 100 1000 500]);

    subplot(1,2,1);
    plot(t, mrcp_pre(cz,:),  'b', 'LineWidth',2, 'DisplayName','PRE'); hold on;
    plot(t, mrcp_post(cz,:), 'r', 'LineWidth',2, 'DisplayName','POST');
    yline(0,'k--','LineWidth',0.8,'HandleVisibility','off');
    xline(0,'k-', 'LineWidth',1.5,'HandleVisibility','off');

    % Shade RP and MP windows
    yl = ylim;
    patch([rp_win(1) rp_win(2) rp_win(2) rp_win(1)],[yl(1) yl(1) yl(2) yl(2)],...
        [0.8 0.9 1.0],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
    patch([mp_win(1) mp_win(2) mp_win(2) mp_win(1)],[yl(1) yl(1) yl(2) yl(2)],...
        [1.0 0.9 0.8],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');

    xlabel('Time re. EMG onset (s)'); ylabel('Amplitude (µV)');
    title(sprintf('MRCP — Cz — %s', strrep(label,'_',' ')));
    legend('Location','best'); grid on;
    xlim([-1.5 1.0]);

    % Metrics table
    subplot(1,2,2);
    axis off;
    txt = {
        sprintf('RP (−0.5→0s):'),
        sprintf('  PRE  = %.3f µV', rp_pre),
        sprintf('  POST = %.3f µV', rp_post),
        sprintf('  Δ    = %+.3f µV', rp_post-rp_pre),
        '',
        sprintf('MP (0→0.1s):'),
        sprintf('  PRE  = %.3f µV', mp_pre),
        sprintf('  POST = %.3f µV', mp_post),
        sprintf('  Δ    = %+.3f µV', mp_post-mp_pre),
    };
    for i = 1:length(txt)
        text(0.05, 0.9 - (i-1)*0.09, txt{i}, 'Units','normalized', ...
            'FontSize',10, 'FontName','Courier');
    end
    title('Metrics at Cz');

    sgtitle(sprintf('MRCP — %s', strrep(label,'_',' ')));
end


function fig = plot_mrcp_grid(mrcp_pre, mrcp_post, t, label, ch_names)
% All 8 channels grid

    fig = figure('Visible','off','Position',[100 100 1400 700]);
    col_pre  = [0.2 0.4 0.8];
    col_post = [0.8 0.2 0.2];

    for ch = 1:8
        subplot(2,4,ch);
        plot(t, mrcp_pre(ch,:),  '-','Color',col_pre, 'LineWidth',1.8,'DisplayName','PRE'); hold on;
        plot(t, mrcp_post(ch,:), '-','Color',col_post,'LineWidth',1.8,'DisplayName','POST');
        yline(0,'k--','LineWidth',0.8,'HandleVisibility','off');
        xline(0,'k-', 'LineWidth',1.2,'HandleVisibility','off');
        title(ch_names{ch},'FontSize',10);
        xlabel('Time (s)'); ylabel('µV');
        xlim([-1.5 1.0]); grid on;
        if ch==1; legend('Location','southwest','FontSize',7); end
    end

    sgtitle(sprintf('MRCP All Channels PRE vs POST — %s', strrep(label,'_',' ')));
end
