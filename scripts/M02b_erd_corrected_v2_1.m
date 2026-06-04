% =========================================================
% MODULE 2b — ERD/ERS, ERP (PRE/POST) & Coherence Analysis
% =========================================================
clear; clc; close all;
eeglab nogui;

ROOT_CLEAN   = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets_clean';
ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';
if ~isfolder(ROOT_RESULTS); mkdir(ROOT_RESULTS); end

FS            = 2400;
TRIAL_SAMPLES = FS * 10;
TRIALS_PER_RUN = 32;

WIN_LEN = round(FS * 0.5);    
OVERLAP = round(FS * 0.45);   
NFFT    = 2048;

BASE_START = 1.5;
BASE_END   = 2.5;

PHASES = struct('prep', [3.5, 5.0], 'exec', [5.0, 7.5], 'recovery', [7.5, 8.5]);
ALPHA = [8,  13];
BETA  = [13, 30];

CH_NAMES = {'C3','C1','Cz','C2','C4','CP1','CP2','Pz'};
MOTOR_CH = [1, 3, 5];   

session_info;
all_results = struct();

for s = 1:size(session_table, 1)
    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);
    folder = fullfile(ROOT_CLEAN, participant, session);
    fprintf('\n=== %s ===\n', label);

    ica_files = dir(fullfile(folder, '*_EEG_ica.set'));
    [~, idx]  = sort({ica_files.name});
    ica_files = ica_files(idx);

    basal_pre_f = ''; basal_post_f = ''; pre_files = {}; post_files = {};
    for f = 1:length(ica_files)
        rn = get_run_number(ica_files(f).name);
        switch rn
            case 1,          basal_pre_f       = ica_files(f).name;
            case {3,4,5,6},  pre_files{end+1}  = ica_files(f).name;
            case {7,8,9,10}, post_files{end+1} = ica_files(f).name;
            case 11,         basal_post_f      = ica_files(f).name;
        end
    end

    % ── ANALYSIS 1: BASAL PSD ──
    if ~isempty(basal_pre_f) && ~isempty(basal_post_f)
        EEG1 = pop_loadset('filename', basal_pre_f,  'filepath', folder);
        EEG2 = pop_loadset('filename', basal_post_f, 'filepath', folder);
        [psd1, ff] = compute_welch_psd(EEG1.data, FS);
        [psd2, ~ ] = compute_welch_psd(EEG2.data, FS);
        
        all_results.(label).psd_pre = psd1; all_results.(label).psd_post = psd2; all_results.(label).freqs = ff;
    end

    % ── ANALYSIS 2: ERD/ERS & COHERENCE ──
    if ~isempty(pre_files) && ~isempty(post_files)
        epochs_pre  = load_epoch_block(pre_files,  folder, TRIAL_SAMPLES);
        epochs_post = load_epoch_block(post_files, folder, TRIAL_SAMPLES);

        % CORREGIDO: ERD real por Trial
        [erd_pre_tf, t_tf, f_tf] = compute_erd_spectrogram_correct(epochs_pre, FS, WIN_LEN, OVERLAP, NFFT, BASE_START, BASE_END, ALPHA, BETA);
        [erd_post_tf, ~, ~]     = compute_erd_spectrogram_correct(epochs_post, FS, WIN_LEN, OVERLAP, NFFT, BASE_START, BASE_END, ALPHA, BETA);
        
        % NUEVO: Coherencia Beta
        [coh_pre, f_coh]  = compute_beta_coherence(epochs_pre, FS, WIN_LEN, OVERLAP, NFFT, BETA);
        [coh_post, ~]     = compute_beta_coherence(epochs_post, FS, WIN_LEN, OVERLAP, NFFT, BETA);

        all_results.(label).erd_pre_tf = erd_pre_tf; all_results.(label).erd_post_tf = erd_post_tf;
        all_results.(label).coherence_pre = coh_pre; all_results.(label).coherence_post = coh_post;
        
        fprintf('    Coherencia Beta C3-Cz: PRE=%.3f | POST=%.3f\n', coh_pre.avg_C3_Cz, coh_post.avg_C3_Cz);
    end

    % ── ANALYSIS 3: ERP PRE Y POST (CORREGIDO) ──
    block_types = {'pre', 'post'};
    files_cell = {pre_files, post_files};
    
    for b = 1:2
        if ~isempty(files_cell{b})
            ep_block = load_epoch_block(files_cell{b}, folder, TRIAL_SAMPLES);
            nyq = FS / 2; [b_lp, a_lp] = butter(4, 20/nyq, 'low');
            epochs_filt = zeros(size(ep_block));
            for tr = 1:size(ep_block,3)
                for ch = 1:8
                    epochs_filt(ch,:,tr) = filtfilt(b_lp, a_lp, ep_block(ch,:,tr));
                end
            end
            stim_s = round(3.0 * FS); pre_win = round(0.5 * FS); post_win = round(2.5 * FS);
            erp_range = (stim_s - pre_win):(stim_s + post_win); t_erp = (-pre_win:post_win) / FS;
            base_idx = t_erp >= -0.5 & t_erp < 0; n_trials = size(epochs_filt, 3);
            erp_bc = zeros(8, length(erp_range), n_trials);
            for tr = 1:n_trials
                seg = epochs_filt(:, erp_range, tr);
                bval = mean(seg(:, base_idx), 2); erp_bc(:,:,tr) = seg - bval;
            end
            all_results.(label).(sprintf('erp_%s', block_types{b})) = mean(erp_bc, 3);
            all_results.(label).t_erp = t_erp;
        end
    end
end
save(fullfile(ROOT_RESULTS, 'all_eeg_results_v4.mat'), 'all_results');


% =========================================================
% HELPER FUNCTIONS
% =========================================================

function run_num = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok); run_num = str2double(tok{1}{1});
    else; run_num = 0; end
end

function [psd, freqs] = compute_welch_psd(data, fs)
    nperseg = fs * 2;
    [~, freqs] = pwelch(data(1,:), nperseg, [], [], fs);
    psd = zeros(size(data,1), length(freqs));
    for ch = 1:size(data,1)
        psd(ch,:) = pwelch(data(ch,:), nperseg, [], [], fs);
    end
end

function bp = band_power(psd, freqs, band)
    idx = freqs >= band(1) & freqs <= band(2);
    bp  = mean(psd(:, idx), 2)';
end

function epochs = load_epoch_block(file_list, folder, trial_samples, trials_per_run)
    data_all = [];
    for f = 1:length(file_list)
        EEG      = pop_loadset('filename', file_list{f}, 'filepath', folder);
        data_all = [data_all, EEG.data];
    end
    n_ch     = size(data_all, 1);
    n_trials = floor(size(data_all,2) / trial_samples);
    epochs   = reshape(data_all(:, 1:n_trials*trial_samples), ...
                       n_ch, trial_samples, n_trials);
end

function [erd_tf, t_tf, f_tf] = compute_erd_spectrogram_correct(epochs, fs, win_len, overlap, nfft, base_start, base_end, alpha, beta)
    n_ch = size(epochs, 1); n_trials = size(epochs, 3);
    [~, f_tf, t_tf, ~] = spectrogram(epochs(1,:,1), win_len, overlap, nfft, fs);
    f_tf = f_tf(:)'; t_tf = t_tf(:)';
    
    idx_alpha = f_tf >= alpha(1) & f_tf <= alpha(2);
    idx_beta  = f_tf >= beta(1)  & f_tf <= beta(2);
    idx_base  = t_tf >= base_start & t_tf <= base_end;
    
    erd_trials_alpha = zeros(n_ch, length(t_tf), n_trials);
    erd_trials_beta  = zeros(n_ch, length(t_tf), n_trials);
    
    for tr = 1:n_trials
        for ch = 1:n_ch
            [~, ~, ~, P] = spectrogram(epochs(ch,:,tr), win_len, overlap, nfft, fs);
            pow_a = mean(P(idx_alpha,:), 1);
            pow_b = mean(P(idx_beta,:), 1);
            
            base_a = mean(pow_a(idx_base));
            base_b = mean(pow_b(idx_base));
            
            % Evitar división por cero
            if base_a == 0; base_a = eps; end
            if base_b == 0; base_b = eps; end
            
            erd_trials_alpha(ch,:,tr) = ((pow_a - base_a) / base_a) * 100;
            erd_trials_beta(ch,:,tr)  = ((pow_b - base_b) / base_b) * 100;
        end
    end
    % Average percents
    erd_tf.alpha = max(min(mean(erd_trials_alpha, 3), 200), -100);
    erd_tf.beta  = max(min(mean(erd_trials_beta, 3), 200), -100);
end

function phases_out = extract_phase_means(erd_tf, t_tf, PHASES)
    ph_names = fieldnames(PHASES);
    bands    = {'alpha','beta'};
    for ph = 1:length(ph_names)
        pname = ph_names{ph};
        pwin  = PHASES.(pname);
        idx   = t_tf >= pwin(1) & t_tf <= pwin(2);
        for ba = 1:length(bands)
            bname = bands{ba};
            phases_out.(pname).(bname) = mean(erd_tf.(bname)(:, idx), 2)';
        end
    end
end

function [coh, f_coh] = compute_beta_coherence(epochs, fs, win_len, overlap, nfft, beta_band)
    % epochs: (canales x muestras x trials)
    n_ch = size(epochs, 1);
    n_trials = size(epochs, 3);
    
    % Obtener frecuencias de referencia usando el primer trial como molde
    [~, f_coh] = mscohere(epochs(1,:,1), epochs(2,:,1), win_len, overlap, nfft, fs);
    idx_beta = f_coh >= beta_band(1) & f_coh <= beta_band(2);
    
    % Inicializar matrices para los pares motoras: 
    % C3 es canal 1, Cz es canal 3, C4 es canal 5
    coh.C3_Cz = zeros(n_trials, 1);
    coh.C4_Cz = zeros(n_trials, 1);
    coh.C3_C4 = zeros(n_trials, 1);
    
    for tr = 1:n_trials
        % Par C3 - Cz
        c_c3_cz = mscohere(epochs(1,:,tr), epochs(3,:,tr), win_len, overlap, nfft, fs);
        coh.C3_Cz(tr) = mean(c_c3_cz(idx_beta));
        
        % Par C4 - Cz
        c_c4_cz = mscohere(epochs(5,:,tr), epochs(3,:,tr), win_len, overlap, nfft, fs);
        coh.C4_Cz(tr) = mean(c_c4_cz(idx_beta));
        
        % Par C3 - C4
        c_c3_c4 = mscohere(epochs(1,:,tr), epochs(5,:,tr), win_len, overlap, nfft, fs);
        coh.C3_C4(tr) = mean(c_c3_c4(idx_beta));
    end
    
    % Promedio final de la coherencia en el bloque (Magnitud al cuadrado de 0 a 1)
    coh.avg_C3_Cz = mean(coh.C3_Cz);
    coh.avg_C4_Cz = mean(coh.C4_Cz);
    coh.avg_C3_C4 = mean(coh.C3_C4);
end

% ── Plot functions ─────────────────────────────────

function fig = plot_erd_timecourse(erd_pre, erd_post, t_tf, label, motor_ch)
% PRE vs POST superimposed — alpha and beta side by side
% Shaded area = difference between PRE and POST (where POST < PRE = more ERD)

    fig = figure('Visible','off','Position',[100 100 1200 500]);

    bands     = {'alpha','beta'};
    band_lbls = {'Alpha (8–13 Hz)','Beta (13–30 Hz)'};
    col_pre   = [0.2 0.4 0.8];
    col_post  = [0.8 0.2 0.2];
    col_diff_more = [0.8 0.2 0.2 0.15];  % POST more negative = more ERD
    col_diff_less = [0.2 0.4 0.8 0.15];  % POST less negative = less ERD

    t_stim   = 3.0;
    t_movin  = 5.0;
    t_movend = 7.5;

    for ba = 1:2
        subplot(1,2,ba);
        bname = bands{ba};

        % Mean across motor channels only
        pre_mean  = mean(erd_pre.(bname)(motor_ch,:),  1);
        post_mean = mean(erd_post.(bname)(motor_ch,:), 1);

        % Shade area between PRE and POST
        % Where POST < PRE (more ERD in POST): red shade
        % Where POST > PRE (less ERD in POST): blue shade
        diff = post_mean - pre_mean;

        % More ERD in POST (diff < 0): shade between curves in red
        mask_more = diff < 0;
        if any(mask_more)
            t_more = t_tf; pre_m = pre_mean; post_m = post_mean;
            pre_m(~mask_more)  = post_mean(~mask_more);  % collapse where not active
            fill([t_more, fliplr(t_more)], ...
                 [pre_m, fliplr(post_m)], ...
                 [0.8 0.2 0.2], 'FaceAlpha', 0.12, 'EdgeColor','none', ...
                 'HandleVisibility','off'); hold on;
        end
        % Less ERD in POST (diff > 0): shade in blue
        mask_less = diff > 0;
        if any(mask_less)
            t_less = t_tf; pre_m = pre_mean; post_m = post_mean;
            pre_m(~mask_less)  = post_mean(~mask_less);
            fill([t_less, fliplr(t_less)], ...
                 [pre_m, fliplr(post_m)], ...
                 [0.2 0.4 0.8], 'FaceAlpha', 0.12, 'EdgeColor','none', ...
                 'HandleVisibility','off'); hold on;
        end

        % Mean lines — on top of shading
        plot(t_tf, pre_mean,  '-', 'Color', col_pre,  'LineWidth', 2.5, ...
            'DisplayName', 'PRE');  hold on;
        plot(t_tf, post_mean, '-', 'Color', col_post, 'LineWidth', 2.5, ...
            'DisplayName', 'POST');

        % Zero line
        yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility','off');

        % Get y limits after plotting
        yl = ylim;

        % Shade execution window (light green)
        patch([t_movin t_movend t_movend t_movin], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.85 0.95 0.85], 'FaceAlpha',0.15,'EdgeColor','none',...
              'HandleVisibility','off');

        % Event markers
        xline(t_stim,   '--', 'Color',[0.6 0.1 0.6], 'LineWidth',1.2, ...
            'Label','Stimulus','LabelVerticalAlignment','bottom',...
            'FontSize',8,'HandleVisibility','off');
        xline(t_movin,  '--', 'Color',[0.1 0.6 0.1], 'LineWidth',1.2, ...
            'Label','MovInit','LabelVerticalAlignment','bottom',...
            'FontSize',8,'HandleVisibility','off');
        xline(t_movend, '--', 'Color',[0.8 0.1 0.1], 'LineWidth',1.2, ...
            'Label','MovEnd','LabelVerticalAlignment','bottom',...
            'FontSize',8,'HandleVisibility','off');

        xlabel('Time (s)');
        ylabel('ERD/ERS (%)');
        title(band_lbls{ba});
        legend('Location','best','FontSize',9);
        xlim([0 10]); grid on;
    end

    sgtitle(sprintf('ERD/ERS Time Course (C3+Cz+C4 mean) — %s', ...
        strrep(label,'_',' ')));
end


function fig = plot_erd_bars(erd_pre, erd_post, label)
    ph_names  = {'prep','exec','recovery'};
    ph_labels = {'Preparation','Execution','Recovery'};
    band_names = {'alpha','beta'};
    band_labels = {'Alpha (8–13 Hz)','Beta (13–30 Hz)'};
    MOTOR_CH = [1, 3, 5];

    fig = figure('Visible','off','Position',[100 100 1000 500]);

    for ba = 1:2
        subplot(1,2,ba);
        bname     = band_names{ba};
        pre_vals  = zeros(1,3);
        post_vals = zeros(1,3);
        for ph = 1:3
            pre_vals(ph)  = mean(erd_pre.(ph_names{ph}).(bname)(MOTOR_CH));
            post_vals(ph) = mean(erd_post.(ph_names{ph}).(bname)(MOTOR_CH));
        end

        x = 1:3;
        b1 = bar(x-0.2, pre_vals,  0.35, 'FaceColor',[0.2 0.4 0.8],'DisplayName','PRE');
        hold on;
        b2 = bar(x+0.2, post_vals, 0.35, 'FaceColor',[0.8 0.2 0.2],'DisplayName','POST');
        yline(0,'k','LineWidth',1,'HandleVisibility','off');

        % Value labels
        for ph = 1:3
            text(ph-0.2, pre_vals(ph)  + sign(pre_vals(ph))*1.5, ...
                sprintf('%.1f',pre_vals(ph)), 'HorizontalAlignment','center','FontSize',8);
            text(ph+0.2, post_vals(ph) + sign(post_vals(ph))*1.5, ...
                sprintf('%.1f',post_vals(ph)),'HorizontalAlignment','center','FontSize',8);
            % Delta annotation
            delta = post_vals(ph) - pre_vals(ph);
            y_ann = max(abs(pre_vals(ph)),abs(post_vals(ph))) + 4;
            if pre_vals(ph) < 0 || post_vals(ph) < 0
                y_ann = -y_ann - 4;
            end
            text(ph, y_ann, sprintf('Δ%+.1f',delta), ...
                'HorizontalAlignment','center','FontSize',7,...
                'Color',[0.3 0.3 0.3],'FontWeight','bold');
        end

        set(gca,'XTick',x,'XTickLabel',ph_labels);
        ylabel('ERD/ERS (%)');
        title(band_labels{ba});
        legend('Location','best');
        grid on;
    end

    sgtitle(sprintf('ERD/ERS PRE vs POST (C3+Cz+C4) — %s', strrep(label,'_',' ')));
end


function fig = plot_erd_channels(erd_pre, erd_post, t_tf, label, ch_names)
% All 8 channels grid — beta band only (most relevant for motor)

    fig = figure('Visible','off','Position',[100 100 1400 700]);

    col_pre  = [0.2 0.4 0.8];
    col_post = [0.8 0.2 0.2];

    for ch = 1:8
        subplot(2,4,ch);
        plot(t_tf, erd_pre.beta(ch,:),  '-', 'Color', col_pre,  ...
            'LineWidth', 1.8, 'DisplayName','PRE');
        hold on;
        plot(t_tf, erd_post.beta(ch,:), '-', 'Color', col_post, ...
            'LineWidth', 1.8, 'DisplayName','POST');
        yline(0,'k--','LineWidth',0.8,'HandleVisibility','off');
        xline(3.0,'--','Color',[0.6 0.1 0.6],'LineWidth',0.8,'HandleVisibility','off');
        xline(5.0,'--','Color',[0.1 0.6 0.1],'LineWidth',0.8,'HandleVisibility','off');
        xline(7.5,'--','Color',[0.8 0.1 0.1],'LineWidth',0.8,'HandleVisibility','off');
        title(ch_names{ch},'FontSize',10);
        xlabel('Time (s)'); ylabel('ERD (%)');
        xlim([0 10]); grid on;
        if ch == 1; legend('Location','southwest','FontSize',7); end
    end

    sgtitle(sprintf('Beta ERD/ERS — All Channels PRE vs POST — %s', ...
        strrep(label,'_',' ')));
end
