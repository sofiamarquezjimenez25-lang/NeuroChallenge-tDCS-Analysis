% =========================================================
% MODULE 6 — Group Statistical Analysis
% NeuroChallenge: EEG + EMG + tDCS
%
% Loads pre-computed results from all modules and applies
% Wilcoxon signed-rank test: delta_tDCS vs delta_sham
%
% Input:  .mat files from results folder
% Output: statistical summary + group figures
% =========================================================

clear; clc; close all;

ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';

% ── Load pre-computed results ──────────────────────────────
eeg  = load(fullfile(ROOT_RESULTS, 'all_eeg_results_v4.mat'));
emg  = load(fullfile(ROOT_RESULTS, 'all_emg_results.mat'));
mrcp = load(fullfile(ROOT_RESULTS, 'all_mrcp.mat'));
cmc  = load(fullfile(ROOT_RESULTS, 'all_cmc.mat'));
lat  = load(fullfile(ROOT_RESULTS, 'all_latency_v2.mat'));

% ── Session table with order info ─────────────────────────
% Format: {participant, session_tDCS, session_sham}
% Add all 8 participants here — update order for P10, P12, P13
% once confirmed
SESSION_ORDER = {
    '5',  '001', '002';   % P5:  S001=tDCS, S002=sham  (UPDATE if wrong)
    '6',  '001', '002';   % P6:  S001=tDCS, S002=sham  (UPDATE if wrong)
    '7',  '001', '002';   % P7:  S001=tDCS, S002=sham
    '8',  '002', '001';   % P8:  S002=tDCS, S001=sham
    '9',  '001', '002';   % P9:  UPDATE once confirmed
    '10', '001', '002';   % P10: UPDATE once confirmed
    '12', '001', '002';   % P12: UPDATE once confirmed
    '13', '001', '002';   % P13: UPDATE once confirmed
};

MOTOR_CH = [1, 3, 5];  % C3, Cz, C4

% ── Extract deltas per participant ─────────────────────────
n_part = size(SESSION_ORDER, 1);

% Pre-allocate — each row = one participant
metrics = struct();
metric_names = {'erd_beta_exec', 'iemg_mean', 'cci', ...
                'erd_alpha_exec', 'erd_beta_prep', 'psd_beta_ratio', ...
                'mp_amplitude', 'cmc_beta', 'delay_cm'};

for m = 1:length(metric_names)
    metrics.(metric_names{m}).tdcs = nan(n_part, 1);
    metrics.(metric_names{m}).sham = nan(n_part, 1);
end

for p = 1:n_part
    pid      = SESSION_ORDER{p, 1};
    sess_tdc = SESSION_ORDER{p, 2};
    sess_sha = SESSION_ORDER{p, 3};

    % Label format used in all result structs
    % e.g. P7_001_stim or P7_002_sham
    % We need to find the correct label regardless of _stim/_sham suffix
    lbl_tdc = find_label(eeg.all_results, pid, sess_tdc);
    lbl_sha = find_label(eeg.all_results, pid, sess_sha);

    if isempty(lbl_tdc) || isempty(lbl_sha)
        fprintf('  [WARNING] Missing data for P%s\n', pid);
        continue;
    end

    % ── ERD beta execution ────────────────────────────────
    % delta = POST - PRE (less negative POST = less ERD = marker of tDCS)
    try
        pre_t  = mean(eeg.all_results.(lbl_tdc).erd_pre_phases.exec.beta(MOTOR_CH));
        post_t = mean(eeg.all_results.(lbl_tdc).erd_post_phases.exec.beta(MOTOR_CH));
        pre_s  = mean(eeg.all_results.(lbl_sha).erd_pre_phases.exec.beta(MOTOR_CH));
        post_s = mean(eeg.all_results.(lbl_sha).erd_post_phases.exec.beta(MOTOR_CH));
        metrics.erd_beta_exec.tdcs(p) = post_t - pre_t;
        metrics.erd_beta_exec.sham(p) = post_s - pre_s;
    catch; end

    % ── ERD alpha execution ───────────────────────────────
    try
        pre_t  = mean(eeg.all_results.(lbl_tdc).erd_pre_phases.exec.alpha(MOTOR_CH));
        post_t = mean(eeg.all_results.(lbl_tdc).erd_post_phases.exec.alpha(MOTOR_CH));
        pre_s  = mean(eeg.all_results.(lbl_sha).erd_pre_phases.exec.alpha(MOTOR_CH));
        post_s = mean(eeg.all_results.(lbl_sha).erd_post_phases.exec.alpha(MOTOR_CH));
        metrics.erd_alpha_exec.tdcs(p) = post_t - pre_t;
        metrics.erd_alpha_exec.sham(p) = post_s - pre_s;
    catch; end

    % ── ERD beta preparation ──────────────────────────────
    try
        pre_t  = mean(eeg.all_results.(lbl_tdc).erd_pre_phases.prep.beta(MOTOR_CH));
        post_t = mean(eeg.all_results.(lbl_tdc).erd_post_phases.prep.beta(MOTOR_CH));
        pre_s  = mean(eeg.all_results.(lbl_sha).erd_pre_phases.prep.beta(MOTOR_CH));
        post_s = mean(eeg.all_results.(lbl_sha).erd_post_phases.prep.beta(MOTOR_CH));
        metrics.erd_beta_prep.tdcs(p) = post_t - pre_t;
        metrics.erd_beta_prep.sham(p) = post_s - pre_s;
    catch; end

    % ── PSD beta ratio (POST/PRE) ─────────────────────────
    try
        ratio_t = mean(eeg.all_results.(lbl_tdc).beta_ratio(MOTOR_CH));
        ratio_s = mean(eeg.all_results.(lbl_sha).beta_ratio(MOTOR_CH));
        metrics.psd_beta_ratio.tdcs(p) = ratio_t;
        metrics.psd_beta_ratio.sham(p) = ratio_s;
    catch; end

    % ── iEMG mean across muscles ──────────────────────────
    try
        pre_t  = mean(mean(emg.all_emg_results.(lbl_tdc).iemg_pre,  2));
        post_t = mean(mean(emg.all_emg_results.(lbl_tdc).iemg_post, 2));
        pre_s  = mean(mean(emg.all_emg_results.(lbl_sha).iemg_pre,  2));
        post_s = mean(mean(emg.all_emg_results.(lbl_sha).iemg_post, 2));
        metrics.iemg_mean.tdcs(p) = post_t - pre_t;
        metrics.iemg_mean.sham(p) = post_s - pre_s;
    catch; end

    % ── CCI ───────────────────────────────────────────────
    try
        pre_t  = mean(emg.all_emg_results.(lbl_tdc).cci_pre);
        post_t = mean(emg.all_emg_results.(lbl_tdc).cci_post);
        pre_s  = mean(emg.all_emg_results.(lbl_sha).cci_pre);
        post_s = mean(emg.all_emg_results.(lbl_sha).cci_post);
        metrics.cci.tdcs(p) = post_t - pre_t;
        metrics.cci.sham(p) = post_s - pre_s;
    catch; end

    % ── MP amplitude MRCP (Cz) ────────────────────────────
    try
        metrics.mp_amplitude.tdcs(p) = ...
            mrcp.all_mrcp.(lbl_tdc).metrics_post.mp_cz - ...
            mrcp.all_mrcp.(lbl_tdc).metrics_pre.mp_cz;
        metrics.mp_amplitude.sham(p) = ...
            mrcp.all_mrcp.(lbl_sha).metrics_post.mp_cz - ...
            mrcp.all_mrcp.(lbl_sha).metrics_pre.mp_cz;
    catch; end

    % ── CMC beta (mean C3+C4, mean EMG channels) ──────────
    try
        b_t = mean(mean(cmc.all_cmc.(lbl_tdc).beta_post - ...
                        cmc.all_cmc.(lbl_tdc).beta_pre, 2));
        b_s = mean(mean(cmc.all_cmc.(lbl_sha).beta_post - ...
                        cmc.all_cmc.(lbl_sha).beta_pre, 2));
        metrics.cmc_beta.tdcs(p) = b_t;
        metrics.cmc_beta.sham(p) = b_s;
    catch; end

    % ── Corticomotor delay ────────────────────────────────
    try
        metrics.delay_cm.tdcs(p) = ...
            lat.all_latency.(lbl_tdc).mean_delay_post - ...
            lat.all_latency.(lbl_tdc).mean_delay_pre;
        metrics.delay_cm.sham(p) = ...
            lat.all_latency.(lbl_sha).mean_delay_post - ...
            lat.all_latency.(lbl_sha).mean_delay_pre;
    catch; end
end

% ── Statistical analysis ───────────────────────────────────
fprintf('\n=== GROUP STATISTICAL ANALYSIS (n=%d) ===\n\n', n_part);

PRIMARY   = {'erd_beta_exec', 'iemg_mean', 'cci'};
SECONDARY = {'erd_alpha_exec', 'erd_beta_prep', 'psd_beta_ratio', ...
             'mp_amplitude'};
EXPLORATORY = {'cmc_beta', 'delay_cm'};

labels_print = struct(...
    'erd_beta_exec',  'ERD beta execution (%)', ...
    'iemg_mean',      'iEMG mean (µV·s)', ...
    'cci',            'CCI (%)', ...
    'erd_alpha_exec', 'ERD alpha execution (%)', ...
    'erd_beta_prep',  'ERD beta preparation (%)', ...
    'psd_beta_ratio', 'PSD beta ratio POST/PRE', ...
    'mp_amplitude',   'MP amplitude Cz (µV)', ...
    'cmc_beta',       'CMC beta delta', ...
    'delay_cm',       'Corticomotor delay delta (s)');

groups = {PRIMARY, SECONDARY, EXPLORATORY};
group_names = {'PRIMARY (Wilcoxon)', 'SECONDARY (effect size)', 'EXPLORATORY'};

stat_results = struct();

for g = 1:3
    fprintf('--- %s ---\n', group_names{g});
    for m = 1:length(groups{g})
        fn = groups{g}{m};
        dt = metrics.(fn).tdcs;
        ds = metrics.(fn).sham;

        % Only use participants with both values
        valid = ~isnan(dt) & ~isnan(ds);
        n_v   = sum(valid);
        dt_v  = dt(valid);
        ds_v  = ds(valid);

        med_t = median(dt_v); iqr_t = iqr(dt_v);
        med_s = median(ds_v); iqr_s = iqr(ds_v);

        if g == 1 && n_v >= 4
            % Wilcoxon signed-rank
            [p_val, ~, stats] = signrank(dt_v, ds_v);
            W = stats.signedrank;
            % Effect size r = Z / sqrt(N)
            % Z approximated from p-value (two-tailed)
            Z = norminv(p_val/2);
            r = abs(Z) / sqrt(n_v);
            sig = '';
            if p_val < 0.05; sig = ' *'; end
            if p_val < 0.01; sig = ' **'; end
            fprintf('  %-30s  tDCS: %.3f [%.3f]  sham: %.3f [%.3f]  W=%g p=%.3f r=%.2f%s  (n=%d)\n', ...
                labels_print.(fn), med_t, iqr_t, med_s, iqr_s, W, p_val, r, sig, n_v);
            stat_results.(fn) = struct('med_tdcs',med_t,'iqr_tdcs',iqr_t,...
                'med_sham',med_s,'iqr_sham',iqr_s,'W',W,'p',p_val,'r',r,'n',n_v);
        else
            % Effect size only (Cohen's d on differences)
            diff = dt_v - ds_v;
            d = mean(diff) / (std(diff) + 1e-10);
            fprintf('  %-30s  tDCS: %.3f [%.3f]  sham: %.3f [%.3f]  d=%.2f  (n=%d)\n', ...
                labels_print.(fn), med_t, iqr_t, med_s, iqr_s, d, n_v);
            stat_results.(fn) = struct('med_tdcs',med_t,'iqr_tdcs',iqr_t,...
                'med_sham',med_s,'iqr_sham',iqr_s,'d',d,'n',n_v);
        end
    end
    fprintf('\n');
end

% ── Figures ───────────────────────────────────────────────

% Figure 1: Raincloud plot for primary metrics
fig1 = plot_raincloud(metrics, PRIMARY, labels_print, n_part);
saveas(fig1, fullfile(ROOT_RESULTS, 'GROUP_primary_raincloud.png'));
close(fig1);

% Figure 2: Heatmap delta tDCS - delta sham across all metrics
fig2 = plot_heatmap(metrics, metric_names, labels_print, SESSION_ORDER, n_part);
saveas(fig2, fullfile(ROOT_RESULTS, 'GROUP_heatmap_deltas.png'));
close(fig2);

% Figure 3: Scatter ERD beta vs CCI (delta level)
fig3 = plot_scatter_eeg_emg(metrics, n_part, SESSION_ORDER);
saveas(fig3, fullfile(ROOT_RESULTS, 'GROUP_ERD_vs_CCI_scatter.png'));
close(fig3);

save(fullfile(ROOT_RESULTS, 'group_stats.mat'), 'metrics', 'stat_results');
fprintf('=== Module 6 complete ===\n');


% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function lbl = find_label(result_struct, pid, sess)
% Finds the field name in result_struct matching participant and session
% e.g. P7_001_stim or P7_002_sham
    lbl = '';
    fields = fieldnames(result_struct);
    for i = 1:length(fields)
        f = fields{i};
        % Match P<pid>_<sess>_
        if contains(f, sprintf('P%s_%s_', pid, sess))
            lbl = f;
            return;
        end
    end
end


function fig = plot_raincloud(metrics, metric_list, labels, n_part)
% Raincloud plot: individual points + jittered strip for tDCS vs sham

    n = length(metric_list);
    fig = figure('Visible','off','Position',[100 100 1200 400]);

    cols = {[0.2 0.4 0.8], [0.8 0.2 0.2]};  % tDCS blue, sham red

    for m = 1:n
        subplot(1, n, m);
        fn = metric_list{m};
        dt = metrics.(fn).tdcs;
        ds = metrics.(fn).sham;
        valid = ~isnan(dt) & ~isnan(ds);
        dt_v = dt(valid); ds_v = ds(valid);

        hold on;
        % Individual points with jitter
        jit = (rand(sum(valid),1)-0.5)*0.15;
        scatter(ones(sum(valid),1) + jit, dt_v, 40, cols{1}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
        scatter(2*ones(sum(valid),1) + jit, ds_v, 40, cols{2}, 'filled', ...
            'MarkerFaceAlpha', 0.6);

        % Lines connecting paired observations
        for i = 1:sum(valid)
            plot([1+jit(i), 2+jit(i)], [dt_v(i), ds_v(i)], ...
                'Color', [0.6 0.6 0.6 0.4], 'LineWidth', 0.8);
        end

        % Median markers
        plot(1, median(dt_v), 's', 'MarkerFaceColor', cols{1}, ...
            'MarkerEdgeColor','k', 'MarkerSize', 10, 'LineWidth', 1.5);
        plot(2, median(ds_v), 's', 'MarkerFaceColor', cols{2}, ...
            'MarkerEdgeColor','k', 'MarkerSize', 10, 'LineWidth', 1.5);

        yline(0, 'k--', 'LineWidth', 1);
        set(gca, 'XTick', [1,2], 'XTickLabel', {'tDCS','sham'}, 'FontSize', 9);
        ylabel('Δ POST−PRE');
        title(labels.(fn), 'FontSize', 8, 'Interpreter', 'none');
        xlim([0.5 2.5]); grid on;
    end

    sgtitle('Primary Metrics: Δ tDCS vs Δ sham (individual participants)');
end


function fig = plot_heatmap(metrics, metric_list, labels, session_order, n_part)
% Heatmap: participants x metrics, color = delta_tDCS - delta_sham

    n_m = length(metric_list);
    data = nan(n_part, n_m);

    for p = 1:n_part
        for m = 1:n_m
            fn = metric_list{m};
            dt = metrics.(fn).tdcs(p);
            ds = metrics.(fn).sham(p);
            if ~isnan(dt) && ~isnan(ds)
                data(p, m) = dt - ds;
            end
        end
    end

    % Normalize each column to [-1, 1] for visualization
    data_norm = data;
    for m = 1:n_m
        col = data(:,m);
        col_valid = col(~isnan(col));
        if ~isempty(col_valid)
            mx = max(abs(col_valid));
            if mx > 0; data_norm(:,m) = col/mx; end
        end
    end

    fig = figure('Visible','off','Position',[100 100 1000 400]);
    imagesc(data_norm, [-1 1]);
    colormap(redblue(256));
    colorbar;

    xlbls = cellfun(@(f) labels.(f), metric_list, 'UniformOutput', false);
    ylbls = cellfun(@(p) sprintf('P%s', p), session_order(:,1), 'UniformOutput', false);

    set(gca, 'XTick', 1:n_m, 'XTickLabel', xlbls, ...
             'XTickLabelRotation', 35, 'FontSize', 8, ...
             'YTick', 1:n_part, 'YTickLabel', ylbls);
    title('Δ tDCS − Δ sham (normalized) — Blue=tDCS effect, Red=sham effect');
    xlabel('Metric'); ylabel('Participant');
end


function fig = plot_scatter_eeg_emg(metrics, n_part, session_order)
% Scatter: Δ ERD beta vs Δ CCI and Δ ERD beta vs Δ iEMG

    fig = figure('Visible','off','Position',[100 100 1000 450]);

    pairs = {{'erd_beta_exec','cci'}, {'erd_beta_exec','iemg_mean'}};
    xlbl  = 'Δ ERD beta execution (%)';
    ylbls = {'Δ CCI (%)', 'Δ iEMG (µV·s)'};

    for sp = 1:2
        subplot(1,2,sp);
        fn_x = pairs{sp}{1};
        fn_y = pairs{sp}{2};

        % tDCS points
        dt_x = metrics.(fn_x).tdcs;
        dt_y = metrics.(fn_y).tdcs;
        ds_x = metrics.(fn_x).sham;
        ds_y = metrics.(fn_y).sham;

        valid_t = ~isnan(dt_x) & ~isnan(dt_y);
        valid_s = ~isnan(ds_x) & ~isnan(ds_y);

        scatter(dt_x(valid_t), dt_y(valid_t), 60, [0.2 0.4 0.8], 'filled', ...
            'DisplayName', 'tDCS'); hold on;
        scatter(ds_x(valid_s), ds_y(valid_s), 60, [0.8 0.2 0.2], 'd', 'filled', ...
            'DisplayName', 'sham');

        % Participant labels
        pids = session_order(:,1);
        for i = 1:n_part
            if valid_t(i)
                text(dt_x(i)+0.5, dt_y(i), sprintf('P%s',pids{i}), 'FontSize',7,'Color',[0.2 0.4 0.8]);
            end
        end

        xline(0,'k--','LineWidth',0.8,'HandleVisibility','off');
        yline(0,'k--','LineWidth',0.8,'HandleVisibility','off');
        xlabel(xlbl); ylabel(ylbls{sp});
        title(sprintf('%s vs %s', 'ΔERD beta', ylbls{sp}(2:end)));
        legend('Location','best','FontSize',8); grid on;
    end

    sgtitle('EEG–EMG relationship: Δ tDCS (circles) vs Δ sham (diamonds)');
end


function c = redblue(m)
% Red-blue diverging colormap centered at zero
    if nargin < 1; m = 64; end
    top    = [1 0 0];
    mid    = [1 1 1];
    bottom = [0 0.4 0.8];
    n = floor(m/2);
    c = [interp1([0 1], [bottom; mid], linspace(0,1,n)); ...
         interp1([0 1], [mid; top],    linspace(0,1,m-n))];
end
