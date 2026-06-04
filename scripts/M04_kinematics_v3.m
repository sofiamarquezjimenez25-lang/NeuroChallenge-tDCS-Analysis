% =========================================================
% MODULE 4 (corrected) — Kinematics from encoder signal
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Since ev_mov_init/ev_mov_end are fixed protocol timestamps
% (not real movement detection), this script detects movement
% onset and offset from the velocity of the encoder signal.
%
% fe and rud are in raw encoder counts — conversion to degrees
% is not applied because the zero position (neutral wrist) is
% unknown per participant. All metrics are computed in counts/s
% and are valid for PRE vs POST comparison within the same session.
%
% Metrics per trial:
%   1. MT  — real movement duration (onset to offset) in seconds
%   2. Peak velocity — max resultant velocity (counts/s)
%   3. Normalized jerk — smoothness metric (lower = smoother)
%
% Movement detection:
%   - Velocity = gradient of encoder signal
%   - Onset: first sample where |velocity| > baseline_mean + 3*SD
%             sustained for at least 50ms
%   - Offset: first sample after onset where |velocity| drops
%             back below threshold and stays there for 100ms
%   - Search window: ev_target_show to ev_target_show + 7s
%
% Input:  .mat files from mat_files folder
% Output: figures and all_kinematics_v3.mat
% =========================================================

clear; clc; close all;

% ── Paths ──────────────────────────────────────────────────
ROOT_MAT     = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\mat_files';
ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';

if ~isfolder(ROOT_RESULTS); mkdir(ROOT_RESULTS); end

% ── Constants ──────────────────────────────────────────────
FS = 2400;

% Movement detection parameters
VEL_THRESHOLD_SD  = 3.0;   % onset threshold: baseline + N*SD
MIN_ONSET_SAMPLES = round(0.05 * FS);   % 50ms sustained above threshold
MIN_OFFSET_SAMPLES= round(0.10 * FS);   % 100ms sustained below threshold

% Baseline window for velocity threshold (seconds within trial)
BASE_WIN = [1.5, 2.5];   % quiet period before stimulus
S_BASE   = round(BASE_WIN(1) * FS) + 1;
E_BASE   = round(BASE_WIN(2) * FS);

% Search window for movement (seconds within trial)
SEARCH_START = round(3.5 * FS) + 1;  % after stimulus
SEARCH_END   = round(9.5 * FS);      % up to near end of trial

% Valid MT range
MT_MIN = 0.1;
MT_MAX = 5.0;

% ── Session table ──────────────────────────────────────────
session_info;

all_kin = struct();

% ── Main loop ──────────────────────────────────────────────
for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);

    mat_folder = fullfile(ROOT_MAT, participant, session);
    fprintf('\n=== %s ===\n', label);

    if ~isfolder(mat_folder)
        fprintf('  [WARNING] Folder not found: %s\n', mat_folder); continue;
    end

    mat_files = dir(fullfile(mat_folder, '*.mat'));
    [~, idx]  = sort({mat_files.name});
    mat_files = mat_files(idx);

    if isempty(mat_files)
        fprintf('  [WARNING] No .mat files\n'); continue;
    end

    % Classify by run
    pre_files = {}; post_files = {};
    for f = 1:length(mat_files)
        rn = get_run_number(mat_files(f).name);
        switch rn
            case {3,4,5,6},  pre_files{end+1}  = mat_files(f).name;
            case {7,8,9,10}, post_files{end+1} = mat_files(f).name;
        end
    end

    fprintf('  PRE: %d runs  |  POST: %d runs\n', ...
        length(pre_files), length(post_files));

    if isempty(pre_files) || isempty(post_files)
        fprintf('  [WARNING] Missing PRE or POST\n'); continue;
    end

    % ── Compute kinematics ────────────────────────────────
    fprintf('  [1] PRE block...\n');
    kin_pre  = compute_kin_block(pre_files,  mat_folder, FS, ...
        S_BASE, E_BASE, SEARCH_START, SEARCH_END, ...
        VEL_THRESHOLD_SD, MIN_ONSET_SAMPLES, MIN_OFFSET_SAMPLES, ...
        MT_MIN, MT_MAX);

    fprintf('  [2] POST block...\n');
    kin_post = compute_kin_block(post_files, mat_folder, FS, ...
        S_BASE, E_BASE, SEARCH_START, SEARCH_END, ...
        VEL_THRESHOLD_SD, MIN_ONSET_SAMPLES, MIN_OFFSET_SAMPLES, ...
        MT_MIN, MT_MAX);

    % ── Print summary ─────────────────────────────────────
    fprintf('\n  === SUMMARY %s ===\n', label);

    metrics = {'mt', 'vpeak', 'jerk_norm'};
    names   = {'MT (s)', 'Peak velocity (counts/s)', 'Norm. jerk'};

    for m = 1:length(metrics)
        fn = metrics{m};
        pv = kin_pre.(fn);  ov = kin_post.(fn);
        mp = mean(pv,'omitnan'); mo = mean(ov,'omitnan');
        sp = std(pv,'omitnan');  so = std(ov,'omitnan');
        np = sum(~isnan(pv));    no = sum(~isnan(ov));
        delta = mo - mp;
        pct   = 100 * delta / (abs(mp) + 1e-10);
        fprintf('  %-28s PRE=%7.3f±%.3f (n=%d)  POST=%7.3f±%.3f (n=%d)  Δ=%+.3f (%+.1f%%)\n', ...
            names{m}, mp, sp, np, mo, so, no, delta, pct);
    end

    fprintf('  Valid onset detections: PRE=%d/%d  POST=%d/%d\n', ...
        sum(~isnan(kin_pre.mt)),  length(kin_pre.mt), ...
        sum(~isnan(kin_post.mt)), length(kin_post.mt));

    % Store
    all_kin.(label).pre  = kin_pre;
    all_kin.(label).post = kin_post;

    % ── Figures ────────────────────────────────────────────

    fig = plot_mt(kin_pre, kin_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_KIN_MT.png', label)));
    close(fig);

    fig = plot_velocity(kin_pre, kin_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_KIN_velocity.png', label)));
    close(fig);

    fig = plot_jerk(kin_pre, kin_post, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_KIN_jerk.png', label)));
    close(fig);

    fig = plot_velocity_profile(kin_pre, kin_post, FS, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_KIN_vel_profile.png', label)));
    close(fig);

    fprintf('  [OK] Figures saved for %s\n', label);
end

% Summary across sessions
fig = plot_summary(all_kin);
saveas(fig, fullfile(ROOT_RESULTS, 'ALL_KIN_summary.png'));
close(fig);

save(fullfile(ROOT_RESULTS, 'all_kinematics_v3.mat'), 'all_kin');
fprintf('\n=== Module 4 complete ===\n');


% =========================================================
% LOCAL FUNCTIONS
% =========================================================

function rn = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok); rn = str2double(tok{1}{1});
    else
        tok = regexp(fname, 'R(\d+)', 'tokens');
        if ~isempty(tok); rn = str2double(tok{1}{1}); else; rn = -1; end
    end
end


function kin = compute_kin_block(file_list, folder, fs, ...
        s_base, e_base, s_search, e_search, ...
        vel_sd, min_onset, min_offset, mt_min, mt_max)
% Computes MT, peak velocity, and jerk from encoder signal
% using velocity-based movement detection

    kin.mt        = [];
    kin.vpeak     = [];
    kin.jerk_norm = [];
    kin.vel_fe    = {};   % store mean velocity profiles for plotting
    kin.vel_rud   = {};

    for f = 1:length(file_list)
        d = load(fullfile(folder, file_list{f}));

        if ~isfield(d,'fe') || ~isfield(d,'rud')
            fprintf('    [WARNING] Missing fe/rud in %s\n', file_list{f});
            continue;
        end

        fe  = double(d.fe(:));
        rud = double(d.rud(:));

        % Determine number of trials from file
        if isfield(d, 'targets')
            n_trials = length(d.targets);
        else
            n_trials = 32;  % default
        end

        trial_samples = floor(length(fe) / n_trials);
        if trial_samples < 1000
            fprintf('    [WARNING] Unexpected signal length in %s\n', file_list{f});
            continue;
        end

        % Epoch into trials
        fe_ep  = reshape(fe(1:n_trials*trial_samples),  trial_samples, n_trials);
        rud_ep = reshape(rud(1:n_trials*trial_samples), trial_samples, n_trials);

        mt_run    = nan(n_trials, 1);
        vpeak_run = nan(n_trials, 1);
        jerk_run  = nan(n_trials, 1);

        for tr = 1:n_trials
            fe_tr  = fe_ep(:, tr);
            rud_tr = rud_ep(:, tr);

            % Velocity (counts/sample → counts/s)
            vel_fe  = gradient(fe_tr)  * fs;
            vel_rud = gradient(rud_tr) * fs;
            vel_mag = sqrt(vel_fe.^2 + vel_rud.^2);

            % Baseline threshold from quiet period
            s_b = min(s_base, trial_samples);
            e_b = min(e_base, trial_samples);
            if e_b <= s_b; continue; end

            base_vel  = vel_mag(s_b:e_b);
            threshold = mean(base_vel) + vel_sd * std(base_vel);

            % Search window
            s_s = min(s_search, trial_samples);
            e_s = min(e_search, trial_samples);
            if e_s <= s_s; continue; end

            search_vel = vel_mag(s_s:e_s);

            % Detect onset: sustained above threshold
            above      = search_vel > threshold;
            onset_idx  = find_sustained(above, min_onset, true);
            if isnan(onset_idx); continue; end

            % Detect offset: sustained below threshold after onset
            post_onset = above(onset_idx:end);
            offset_rel = find_sustained(post_onset, min_offset, false);
            if isnan(offset_rel)
                offset_idx = length(search_vel);  % use end of search window
            else
                offset_idx = onset_idx + offset_rel - 1;
            end

            % Convert to absolute sample indices
            abs_onset  = s_s + onset_idx  - 1;
            abs_offset = s_s + offset_idx - 1;
            abs_offset = min(abs_offset, trial_samples);

            % MT in seconds
            mt_val = (abs_offset - abs_onset) / fs;
            if mt_val < mt_min || mt_val > mt_max; continue; end
            mt_run(tr) = mt_val;

            % Movement segment
            fe_seg  = fe_tr(abs_onset:abs_offset);
            rud_seg = rud_tr(abs_onset:abs_offset);
            T = mt_val;

            % Peak velocity
            vel_seg = vel_mag(abs_onset:abs_offset);
            vpeak_run(tr) = max(vel_seg);

            % Normalized jerk
            if T > 0.1 && length(fe_seg) > 20
                v_fe  = gradient(fe_seg)  * fs;
                v_rud = gradient(rud_seg) * fs;
                a_fe  = gradient(v_fe)    * fs;
                a_rud = gradient(v_rud)   * fs;
                j_fe  = gradient(a_fe)    * fs;
                j_rud = gradient(a_rud)   * fs;

                jrk_int = trapz(j_fe.^2 + j_rud.^2) / fs;

                D_fe  = fe_seg(end)  - fe_seg(1);
                D_rud = rud_seg(end) - rud_seg(1);
                D = sqrt(D_fe^2 + D_rud^2);

                if D > 0.5
                    jerk_run(tr) = sqrt(0.5 * jrk_int) * (T^2.5 / D);
                end
            end
        end

        kin.mt        = [kin.mt;        mt_run];
        kin.vpeak     = [kin.vpeak;     vpeak_run];
        kin.jerk_norm = [kin.jerk_norm; jerk_run];
    end
end


function idx = find_sustained(binary_vec, min_samples, find_true)
% Finds first index where binary_vec is consistently true (or false)
% for at least min_samples consecutive samples
% find_true=true: find first sustained TRUE
% find_true=false: find first sustained FALSE

    idx = NaN;
    n   = length(binary_vec);

    if find_true
        target = true;
    else
        target = false;
        binary_vec = ~binary_vec;  % invert: find first sustained false = first sustained true of inverted
    end

    count = 0;
    for i = 1:n
        if binary_vec(i) == target
            count = count + 1;
            if count >= min_samples
                idx = i - min_samples + 1;
                return;
            end
        else
            count = 0;
        end
    end
end


% ── Plot functions ──────────────────────────────────────────

function fig = plot_mt(pre, post, label)
    fig = figure('Visible','off','Position',[100 100 1000 420]);

    vp = pre.mt(~isnan(pre.mt));
    vo = post.mt(~isnan(post.mt));

    subplot(1,2,1);
    if ~isempty(vp) || ~isempty(vo)
        all_v = [vp;vo];
        edges = linspace(max(0,min(all_v)*0.8), max(all_v)*1.1, 25);
        histogram(vp, edges,'FaceColor',[0.2 0.4 0.8],'FaceAlpha',0.6,...
            'DisplayName',sprintf('PRE (n=%d)',length(vp)));
        hold on;
        histogram(vo, edges,'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.6,...
            'DisplayName',sprintf('POST (n=%d)',length(vo)));
        if ~isempty(vp); xline(mean(vp),'b--','LineWidth',2,'HandleVisibility','off'); end
        if ~isempty(vo); xline(mean(vo),'r--','LineWidth',2,'HandleVisibility','off'); end
        xlabel('Movement Time (s)'); ylabel('Trials');
        title('MT distribution'); legend('Location','best'); grid on;
    end

    subplot(1,2,2);
    means = [mean(vp,'omitnan'), mean(vo,'omitnan')];
    stds  = [std(vp,'omitnan'),  std(vo,'omitnan')];
    if ~any(isnan(means))
        b = bar([1,2],means,0.5); b.FaceColor='flat';
        b.CData = [0.2 0.4 0.8; 0.8 0.2 0.2];
        hold on;
        errorbar([1,2],means,stds,'k','LineStyle','none','LineWidth',1.5,'CapSize',8);
        delta = means(2)-means(1);
        text(1.5,max(means+stds)*1.12,sprintf('Δ=%+.3fs\n(%+.1f%%)',delta,100*delta/means(1)),...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    end
    set(gca,'XTick',[1,2],'XTickLabel',{'PRE','POST'});
    ylabel('Mean MT (s)'); title('Mean ± SD'); grid on;

    sgtitle(sprintf('Movement Time — %s', strrep(label,'_',' ')));
end


function fig = plot_velocity(pre, post, label)
    fig = figure('Visible','off','Position',[100 100 1000 420]);

    vp = pre.vpeak(~isnan(pre.vpeak));
    vo = post.vpeak(~isnan(post.vpeak));

    subplot(1,2,1);
    if ~isempty(vp) || ~isempty(vo)
        all_v = [vp;vo];
        edges = linspace(0, prctile(all_v,97)*1.1, 25);
        histogram(vp, edges,'FaceColor',[0.2 0.4 0.8],'FaceAlpha',0.6,...
            'DisplayName',sprintf('PRE (n=%d)',length(vp)));
        hold on;
        histogram(vo, edges,'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.6,...
            'DisplayName',sprintf('POST (n=%d)',length(vo)));
        if ~isempty(vp); xline(mean(vp),'b--','LineWidth',2,'HandleVisibility','off'); end
        if ~isempty(vo); xline(mean(vo),'r--','LineWidth',2,'HandleVisibility','off'); end
        xlabel('Peak velocity (counts/s)'); ylabel('Trials');
        title('Peak velocity distribution'); legend('Location','best'); grid on;
    end

    subplot(1,2,2);
    means = [mean(vp,'omitnan'), mean(vo,'omitnan')];
    stds  = [std(vp,'omitnan'),  std(vo,'omitnan')];
    if ~any(isnan(means))
        b = bar([1,2],means,0.5); b.FaceColor='flat';
        b.CData = [0.2 0.4 0.8; 0.8 0.2 0.2];
        hold on;
        errorbar([1,2],means,stds,'k','LineStyle','none','LineWidth',1.5,'CapSize',8);
        delta = means(2)-means(1);
        text(1.5,max(means+stds)*1.12,sprintf('Δ=%+.1f\n(%+.1f%%)',delta,100*delta/means(1)),...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    end
    set(gca,'XTick',[1,2],'XTickLabel',{'PRE','POST'});
    ylabel('Mean peak velocity (counts/s)'); title('Mean ± SD'); grid on;

    sgtitle(sprintf('Peak Velocity — %s', strrep(label,'_',' ')));
end


function fig = plot_jerk(pre, post, label)
    fig = figure('Visible','off','Position',[100 100 1000 420]);

    vp = pre.jerk_norm(~isnan(pre.jerk_norm));
    vo = post.jerk_norm(~isnan(post.jerk_norm));

    subplot(1,2,1);
    if ~isempty(vp) || ~isempty(vo)
        all_v = [vp;vo];
        edges = linspace(0, prctile(all_v,95)*1.1, 25);
        histogram(vp, edges,'FaceColor',[0.2 0.4 0.8],'FaceAlpha',0.6,...
            'DisplayName',sprintf('PRE (n=%d)',length(vp)));
        hold on;
        histogram(vo, edges,'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.6,...
            'DisplayName',sprintf('POST (n=%d)',length(vo)));
        if ~isempty(vp); xline(mean(vp),'b--','LineWidth',2,'HandleVisibility','off'); end
        if ~isempty(vo); xline(mean(vo),'r--','LineWidth',2,'HandleVisibility','off'); end
        xlabel('Normalized jerk'); ylabel('Trials');
        title('Jerk distribution (lower = smoother)');
        legend('Location','best'); grid on;
    end

    subplot(1,2,2);
    means = [mean(vp,'omitnan'), mean(vo,'omitnan')];
    stds  = [std(vp,'omitnan'),  std(vo,'omitnan')];
    if ~any(isnan(means))
        b = bar([1,2],means,0.5); b.FaceColor='flat';
        b.CData = [0.2 0.4 0.8; 0.8 0.2 0.2];
        hold on;
        errorbar([1,2],means,stds,'k','LineStyle','none','LineWidth',1.5,'CapSize',8);
        delta = means(2)-means(1);
        text(1.5,max(means+stds)*1.12,sprintf('Δ=%+.3f\n(%+.1f%%)',delta,100*delta/(means(1)+1e-10)),...
            'HorizontalAlignment','center','FontWeight','bold','FontSize',10);
    end
    set(gca,'XTick',[1,2],'XTickLabel',{'PRE','POST'});
    ylabel('Mean norm. jerk'); title('Mean ± SD'); grid on;

    sgtitle(sprintf('Smoothness — %s', strrep(label,'_',' ')));
end


function fig = plot_velocity_profile(pre, post, fs, label)
% Average velocity profile aligned to movement onset

    fig = figure('Visible','off','Position',[100 100 900 420]);

    % Reconstruct from stored trials if available
    % Use vpeak distribution as proxy for profile shape
    % Since we don't store full segments, show histogram comparison

    % Time axis for a typical movement (0 to 3s)
    % This panel shows the mean velocity magnitude across trials
    % aligned approximately to movement period

    vp = pre.vpeak(~isnan(pre.vpeak));
    vo = post.vpeak(~isnan(post.vpeak));
    mp = pre.mt(~isnan(pre.mt));
    mo = post.mt(~isnan(post.mt));

    subplot(1,2,1);
    % MT vs peak velocity scatter
    valid_p = ~isnan(pre.mt) & ~isnan(pre.vpeak);
    valid_o = ~isnan(post.mt) & ~isnan(post.vpeak);
    scatter(pre.mt(valid_p), pre.vpeak(valid_p), 15, [0.2 0.4 0.8], 'filled', ...
        'MarkerFaceAlpha',0.4,'DisplayName','PRE');
    hold on;
    scatter(post.mt(valid_o), post.vpeak(valid_o), 15, [0.8 0.2 0.2], 'filled', ...
        'MarkerFaceAlpha',0.4,'DisplayName','POST');
    xlabel('MT (s)'); ylabel('Peak velocity (counts/s)');
    title('MT vs Peak Velocity per trial');
    legend('Location','best'); grid on;

    subplot(1,2,2);
    % Boxplot MT and velocity
    data_mt = {mp, mo};
    cols    = {[0.2 0.4 0.8], [0.8 0.2 0.2]};
    lbls    = {'PRE MT','POST MT'};
    for i = 1:2
        v = data_mt{i};
        if isempty(v); continue; end
        q = quantile(v,[0.25 0.5 0.75]);
        iqrv = q(3)-q(1);
        wlo = max(min(v), q(1)-1.5*iqrv);
        whi = min(max(v), q(3)+1.5*iqrv);
        rectangle('Position',[i-0.3,q(1),0.6,q(3)-q(1)],...
            'FaceColor',[cols{i},0.4],'EdgeColor',cols{i},'LineWidth',1.5);
        line([i-0.3,i+0.3],[q(2),q(2)],'Color',cols{i},'LineWidth',2.5);
        line([i,i],[wlo,q(1)],'Color',cols{i},'LineWidth',1.5);
        line([i,i],[q(3),whi],'Color',cols{i},'LineWidth',1.5);
        plot(i,mean(v),'o','MarkerFaceColor',cols{i},'MarkerEdgeColor','k','MarkerSize',8);
    end
    set(gca,'XTick',1:2,'XTickLabel',lbls);
    ylabel('MT (s)'); title('MT boxplot PRE vs POST'); grid on;

    sgtitle(sprintf('Movement Profile — %s', strrep(label,'_',' ')));
end


function fig = plot_summary(all_kin)
    fields = fieldnames(all_kin);
    n = length(fields);

    metrics = {'mt','vpeak','jerk_norm'};
    names   = {'MT (s)','Peak velocity (counts/s)','Norm. jerk'};

    fig = figure('Visible','off','Position',[100 100 1400 500]);

    for m = 1:3
        subplot(1,3,m);
        fn = metrics{m};

        pre_m = nan(n,1); post_m = nan(n,1);
        pre_s = nan(n,1); post_s = nan(n,1);
        lbls  = cell(n,1);

        for i = 1:n
            if isfield(all_kin.(fields{i}).pre, fn)
                pv = all_kin.(fields{i}).pre.(fn);
                ov = all_kin.(fields{i}).post.(fn);
                pre_m(i)  = mean(pv,'omitnan');
                post_m(i) = mean(ov,'omitnan');
                pre_s(i)  = std(pv,'omitnan');
                post_s(i) = std(ov,'omitnan');
            end
            lbls{i} = strrep(fields{i},'_',' ');
        end

        x = 1:n;
        bar(x-0.2, pre_m,  0.35,'FaceColor',[0.2 0.4 0.8],'DisplayName','PRE');
        hold on;
        bar(x+0.2, post_m, 0.35,'FaceColor',[0.8 0.2 0.2],'DisplayName','POST');
        errorbar(x-0.2,pre_m, pre_s, 'k','LineStyle','none','LineWidth',1,'CapSize',5);
        errorbar(x+0.2,post_m,post_s,'k','LineStyle','none','LineWidth',1,'CapSize',5);

        for i = 1:n
            if ~isnan(pre_m(i)) && ~isnan(post_m(i))
                delta = post_m(i)-pre_m(i);
                y_top = max(pre_m(i)+pre_s(i),post_m(i)+post_s(i)) + ...
                        abs(max(pre_m(i),post_m(i)))*0.05;
                text(i,y_top,sprintf('%+.2f',delta),...
                    'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
            end
        end

        set(gca,'XTick',x,'XTickLabel',lbls,'FontSize',7);
        ylabel(names{m}); title(names{m});
        if m==1; legend('Location','best'); end
        grid on;
    end

    sgtitle('Kinematics Summary — All Sessions');
end
