% =========================================================
% MODULE 5a — Corticomuscular Coherence (CMC)
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Computes spectral coherence between EEG motor channels
% (C3, C4) and rectified bipolar EMG during the execution
% window (5.0–7.5 s) across all trials in each block.
%
% CMC = |Cxy(f)|^2 = |Pxy(f)|^2 / (Pxx(f) * Pyy(f))
%
% where Pxy is the cross-spectrum, Pxx and Pyy are the
% auto-spectra of EEG and EMG respectively.
%
% Expected finding: if tDCS increased M1 excitability,
% beta-band CMC (13–30 Hz) should be higher POST vs PRE
% in the tDCS session, with a larger delta than sham.
%
% Method: Welch-based coherence (mscohere) averaged across
% trials, separately for each EEG-EMG channel pair.
%
% Input:
%   - *_EEG_ica.set (sets_clean folder)
%   - *_EMG.set (sets folder)
%
% Output:
%   - CMC spectra figures PRE vs POST per session
%   - Bar chart: mean beta CMC PRE vs POST
%   - Summary metrics saved to all_cmc.mat
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
TRIAL_SAMPLES = FS * 10;       % 24000 samples per 10s trial
TRIALS_PER_RUN = 32;

% Execution window (samples within each trial)
S_EXEC = round(5.0 * FS) + 1;  % sample 12001
E_EXEC = round(7.5 * FS);       % sample 18000
EXEC_LEN = E_EXEC - S_EXEC + 1; % 6000 samples = 2.5 s

% Coherence estimation parameters
% Window: 0.5 s, 50% overlap → good freq resolution at 2 Hz
WIN_COH  = round(0.5 * FS);    % 1200 samples
OVERLAP  = round(0.25 * FS);   % 50% overlap
NFFT_COH = 2048;               % freq resolution ~1.17 Hz

% Frequency bands
BETA = [13, 30];

% EEG channels of interest
CH_NAMES = {'C3','C1','Cz','C2','C4','CP1','CP2','Pz'};
CH_C3 = 1;
CH_C4 = 5;
EEG_CMC_CH = [CH_C3, CH_C4];   % C3 and C4

% EMG muscles
MUSCLE_NAMES = {'M1','M2','M3','M4'};

% ── Session table ──────────────────────────────────────────
session_info;

all_cmc = struct();

% ── Main loop ──────────────────────────────────────────────
for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);

    eeg_folder = fullfile(ROOT_CLEAN, participant, session);
    emg_folder = fullfile(ROOT_SETS,  participant, session);

    fprintf('\n=== %s ===\n', label);

    % Classify files
    ica_files     = dir(fullfile(eeg_folder, '*_EEG_ica.set'));
    emg_files_all = dir(fullfile(emg_folder, '*_EMG.set'));
    [~, idx] = sort({ica_files.name});     ica_files     = ica_files(idx);
    [~, idx] = sort({emg_files_all.name}); emg_files_all = emg_files_all(idx);

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
        fprintf('  [WARNING] Missing files, skipping\n');
        continue;
    end

    % ── Load and preprocess EEG blocks ────────────────────
    fprintf('  Loading EEG...\n');
    [eeg_pre_ep,  n_eeg_pre]  = load_eeg_epochs(pre_eeg,  eeg_folder, TRIAL_SAMPLES, S_EXEC, E_EXEC);
    [eeg_post_ep, n_eeg_post] = load_eeg_epochs(post_eeg, eeg_folder, TRIAL_SAMPLES, S_EXEC, E_EXEC);

    % ── Load and preprocess EMG blocks ────────────────────
    fprintf('  Loading EMG...\n');
    [emg_pre_ep,  n_emg_pre]  = load_emg_epochs(pre_emg,  emg_folder, FS, TRIAL_SAMPLES, S_EXEC, E_EXEC);
    [emg_post_ep, n_emg_post] = load_emg_epochs(post_emg, emg_folder, FS, TRIAL_SAMPLES, S_EXEC, E_EXEC);

    n_pre  = min(n_eeg_pre,  n_emg_pre);
    n_post = min(n_eeg_post, n_emg_post);

    fprintf('  Matched trials: PRE=%d  POST=%d\n', n_pre, n_post);

    % ── Compute CMC ───────────────────────────────────────
    fprintf('  Computing CMC...\n');

    % CMC structure: (eeg_ch x emg_ch x freq_bins)
    [cmc_pre,  freqs] = compute_cmc_block(eeg_pre_ep(:,:,1:n_pre),   emg_pre_ep(:,:,1:n_pre), ...
        FS, WIN_COH, OVERLAP, NFFT_COH, EEG_CMC_CH);
    [cmc_post, ~    ] = compute_cmc_block(eeg_post_ep(:,:,1:n_post), emg_post_ep(:,:,1:n_post), ...
        FS, WIN_COH, OVERLAP, NFFT_COH, EEG_CMC_CH);

    % Extract beta-band mean CMC
    idx_beta = freqs >= BETA(1) & freqs <= BETA(2);

    beta_pre  = squeeze(mean(cmc_pre(:,:,idx_beta),  3));  % (eeg_ch x emg_ch)
    beta_post = squeeze(mean(cmc_post(:,:,idx_beta), 3));

    % Print summary
    eeg_lbls = {'C3','C4'};
    fprintf('  Beta CMC summary (mean across EMG channels):\n');
    for ec = 1:length(EEG_CMC_CH)
        pre_mean  = mean(beta_pre(ec,:));
        post_mean = mean(beta_post(ec,:));
        fprintf('    %s: PRE=%.4f  POST=%.4f  delta=%+.4f\n', ...
            eeg_lbls{ec}, pre_mean, post_mean, post_mean - pre_mean);
    end

    % Store results
    all_cmc.(label).cmc_pre   = cmc_pre;
    all_cmc.(label).cmc_post  = cmc_post;
    all_cmc.(label).freqs     = freqs;
    all_cmc.(label).beta_pre  = beta_pre;
    all_cmc.(label).beta_post = beta_post;
    all_cmc.(label).n_pre     = n_pre;
    all_cmc.(label).n_post    = n_post;

    % ── Figures ───────────────────────────────────────────

    % Figure 1: CMC spectrum PRE vs POST (mean across EMG channels)
    fig = plot_cmc_spectrum(cmc_pre, cmc_post, freqs, BETA, label, eeg_lbls, MUSCLE_NAMES, n_pre, n_post);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_CMC_spectrum.png', label)));
    close(fig);

    % Figure 2: Beta CMC bar chart PRE vs POST
    fig = plot_beta_bars(beta_pre, beta_post, label, eeg_lbls, MUSCLE_NAMES);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_CMC_beta_bars.png', label)));
    close(fig);

    fprintf('  [OK] CMC figures saved for %s\n', label);
end

save(fullfile(ROOT_RESULTS, 'all_cmc.mat'), 'all_cmc');
fprintf('\n=== Module 5a complete ===\n');


% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function run_num = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok)
        run_num = str2double(tok{1}{1});
    else
        tok = regexp(fname, 'run(\d+)', 'tokens', 'ignorecase');
        if ~isempty(tok); run_num = str2double(tok{1}{1});
        else; run_num = 0; end
    end
end


function [epochs, n_trials] = load_eeg_epochs(file_list, folder, trial_samples, s_exec, e_exec)
% Loads EEG runs, segments into trials, extracts execution window
% Output: epochs (8 x exec_samples x n_trials)

    exec_len = e_exec - s_exec + 1;
    data_all = [];
    for f = 1:length(file_list)
        EEG      = pop_loadset('filename', file_list{f}, 'filepath', folder);
        data_all = [data_all, EEG.data]; %#ok<AGROW>
    end
    n_ch     = size(data_all, 1);
    n_trials = floor(size(data_all, 2) / trial_samples);
    data_ep  = reshape(data_all(:, 1:n_trials*trial_samples), n_ch, trial_samples, n_trials);
    epochs   = data_ep(:, s_exec:e_exec, :);  % (ch x exec_samples x trials)
    n_trials = size(epochs, 3);
end


function [epochs, n_trials] = load_emg_epochs(file_list, folder, fs, trial_samples, s_exec, e_exec)
% Loads EMG runs, bipolarizes, bandpass, rectifies, extracts execution window
% Output: epochs (4 x exec_samples x n_trials) — rectified bipolar EMG

    nyq = fs / 2;
    [b_bp, a_bp] = butter(4, [20 500]/nyq, 'bandpass');

    data_all = [];
    for f = 1:length(file_list)
        EEG      = pop_loadset('filename', file_list{f}, 'filepath', folder);
        raw      = EEG.data;

        % Bipolarize
        bipolar = zeros(4, size(raw, 2));
        for m = 1:4
            bipolar(m,:) = raw(2*m-1,:) - raw(2*m,:);
        end

        % Bandpass
        for m = 1:4
            bipolar(m,:) = filtfilt(b_bp, a_bp, bipolar(m,:));
        end

        % Rectify (for CMC use rectified EMG, not envelope)
        bipolar = abs(bipolar);

        data_all = [data_all, bipolar]; %#ok<AGROW>
    end

    n_trials = floor(size(data_all, 2) / trial_samples);
    data_ep  = reshape(data_all(:, 1:n_trials*trial_samples), 4, trial_samples, n_trials);
    epochs   = data_ep(:, s_exec:e_exec, :);
    n_trials = size(epochs, 3);
end


function [cmc, freqs] = compute_cmc_block(eeg_ep, emg_ep, fs, win_len, overlap, nfft, eeg_ch_idx)
% Correct CMC: average coherences per trial (not spectra)

    n_eeg    = length(eeg_ch_idx);
    n_emg    = size(emg_ep, 1);
    n_trials = size(eeg_ep, 3);

    % Frequencies from first call
    [~, freqs] = mscohere(eeg_ep(eeg_ch_idx(1),:,1), emg_ep(1,:,1), ...
        hann(win_len), overlap, nfft, fs);
    n_freqs = length(freqs);

    % Accumulate COHERENCE (not spectra) per trial
    cmc_sum = zeros(n_eeg, n_emg, n_freqs);

    for tr = 1:n_trials
        for ec = 1:n_eeg
            eeg_sig = eeg_ep(eeg_ch_idx(ec), :, tr);
            for em = 1:n_emg
                emg_sig = emg_ep(em, :, tr);
                % mscohere directly returns true coherence [0,1]
                coh = mscohere(eeg_sig, emg_sig, ...
                    hann(win_len), overlap, nfft, fs);
                cmc_sum(ec, em, :) = squeeze(cmc_sum(ec,em,:))' + coh';
            end
        end
    end

    % Average consistency between trials
    cmc = cmc_sum / n_trials;
end


% ── Plot functions ──────────────────────────────────────────

function fig = plot_cmc_spectrum(cmc_pre, cmc_post, freqs, beta, label, eeg_lbls, muscle_names, n_pre, n_post)
% CMC spectrum PRE vs POST, one subplot per EEG channel
% Each curve is the mean CMC across the 4 EMG channels

    fig = figure('Visible','off','Position',[100 100 1000 450]);

    n_eeg = size(cmc_pre, 1);
    cols_pre  = [0.2 0.4 0.8];
    cols_post = [0.8 0.2 0.2];

    for ec = 1:n_eeg
        subplot(1, n_eeg, ec);

        % Mean across EMG channels
        cmc_pre_mean  = squeeze(mean(cmc_pre(ec,:,:),  2));
        cmc_post_mean = squeeze(mean(cmc_post(ec,:,:), 2));

        % Individual EMG channels (thin)
        for em = 1:size(cmc_pre, 2)
            plot(freqs, squeeze(cmc_pre(ec,em,:)),  '-', 'Color', [cols_pre,  0.2], 'LineWidth', 0.8, 'HandleVisibility','off'); hold on;
            plot(freqs, squeeze(cmc_post(ec,em,:)), '-', 'Color', [cols_post, 0.2], 'LineWidth', 0.8, 'HandleVisibility','off');
        end

        % Mean (thick)
        plot(freqs, cmc_pre_mean,  '-', 'Color', cols_pre,  'LineWidth', 2.5, 'DisplayName', sprintf('PRE (n=%d)', n_pre));
        plot(freqs, cmc_post_mean, '-', 'Color', cols_post, 'LineWidth', 2.5, 'DisplayName', sprintf('POST (n=%d)', n_post));

        % Beta shading
        ylims = [0, max([cmc_pre_mean(:); cmc_post_mean(:)])*1.3 + 0.001];
        patch([beta(1) beta(2) beta(2) beta(1)], [0 0 ylims(2) ylims(2)], ...
            [0.9 0.95 1.0], 'FaceAlpha', 0.2, 'EdgeColor','none');
        text(mean(beta), ylims(2)*0.95, 'beta', 'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.8]);

        xlim([5 50]);
        ylim(ylims);
        xlabel('Frequency (Hz)');
        ylabel('Coherence');
        title(sprintf('EEG %s — mean EMG M1–M4', eeg_lbls{ec}));
        legend('Location','northeast','FontSize',8);
        grid on;
    end

    sgtitle(sprintf('CMC Spectrum PRE vs POST — %s', strrep(label,'_',' ')));
end


function fig = plot_beta_bars(beta_pre, beta_post, label, eeg_lbls, muscle_names)
% Bar chart: beta CMC per EEG-EMG pair, PRE vs POST

    n_eeg = size(beta_pre, 1);
    n_emg = size(beta_pre, 2);

    fig = figure('Visible','off','Position',[100 100 1000 400]);

    for ec = 1:n_eeg
        subplot(1, n_eeg, ec);

        pre_vals  = beta_pre(ec,:);
        post_vals = beta_post(ec,:);

        x = 1:n_emg;
        b1 = bar(x-0.2, pre_vals,  0.35);
        b1.FaceColor = [0.2 0.4 0.8];
        hold on;
        b2 = bar(x+0.2, post_vals, 0.35);
        b2.FaceColor = [0.8 0.2 0.2];

        % Delta annotations
        for em = 1:n_emg
            delta = post_vals(em) - pre_vals(em);
            y_top = max(pre_vals(em), post_vals(em)) * 1.15 + 0.0005;
            text(em, y_top, sprintf('%+.4f', delta), ...
                'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
        end

        set(gca,'XTick',x,'XTickLabel',muscle_names);
        ylabel('Mean beta coherence');
        title(sprintf('EEG %s vs EMG muscles', eeg_lbls{ec}));
        legend('PRE','POST','Location','best');
        grid on;
    end

    sgtitle(sprintf('Beta CMC (13–30 Hz) — %s', strrep(label,'_',' ')));
end
