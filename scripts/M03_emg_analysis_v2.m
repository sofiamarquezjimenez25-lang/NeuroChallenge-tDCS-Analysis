% =========================================================
% MODULE 3 — EMG Analysis
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Steps:
%   1. Bipolarization (M1=M1L-M1R, etc.)
%   2. Bandpass filter 20-500 Hz
%   3. Rectification + envelope (LP 8 Hz)
%   4. RMS per phase
%   5. iEMG (integrated EMG)
%   6. CCI (co-contraction index, symmetric pairs)
%
% RT removed — use M02c for reaction time (more robust)
% =========================================================

clear; clc; close all;
eeglab nogui;

ROOT_SETS    = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets';
ROOT_RESULTS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\results';

if ~isfolder(ROOT_RESULTS); mkdir(ROOT_RESULTS); end

FS             = 2400;
TRIAL_SAMPLES  = FS * 10;
TRIALS_PER_RUN = 32;

WIN_BASE = [1.5, 2.5];
WIN_PREP = [3.5, 5.0];
WIN_EXEC = [5.0, 7.5];
WIN_REC  = [7.5, 8.5];

EMG_BIPO = {'M1','M2','M3','M4'};

session_info;
all_emg_results = struct();

for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};
    label       = sprintf('P%s_%s_%s', participant, session, sess_type);
    folder      = fullfile(ROOT_SETS, participant, session);

    fprintf('\n=== %s ===\n', label);

    emg_files = dir(fullfile(folder, '*_EMG.set'));
    [~, idx]  = sort({emg_files.name});
    emg_files = emg_files(idx);

    pre_files = {}; post_files = {}; stim_files = {};
    for f = 1:length(emg_files)
        rn = get_run_number(emg_files(f).name);
        switch rn
            case {3,4,5,6},  pre_files{end+1}  = emg_files(f).name;
            case 2,          stim_files{end+1} = emg_files(f).name;
            case {7,8,9,10}, post_files{end+1} = emg_files(f).name;
        end
    end

    if isempty(pre_files) || isempty(post_files)
        fprintf('  [WARNING] Missing PRE or POST\n'); continue;
    end

    fprintf('  PRE: %d runs | POST: %d runs | STIM: %d runs\n', ...
        length(pre_files), length(post_files), length(stim_files));

    env_pre  = process_emg_block(pre_files,  folder, FS, TRIAL_SAMPLES);
    env_post = process_emg_block(post_files, folder, FS, TRIAL_SAMPLES);
    env_stim = [];
    if ~isempty(stim_files)
        env_stim = process_emg_block(stim_files, folder, FS, TRIAL_SAMPLES);
    end

    % ── RMS ───────────────────────────────────────────────
    fprintf('  [1] RMS...\n');
    rms_pre  = compute_rms_phases(env_pre,  FS, WIN_BASE, WIN_PREP, WIN_EXEC, WIN_REC);
    rms_post = compute_rms_phases(env_post, FS, WIN_BASE, WIN_PREP, WIN_EXEC, WIN_REC);

    max_rms  = max([max(rms_pre.exec(:)), max(rms_post.exec(:))]) + 1e-10;
    rms_pre  = normalize_rms(rms_pre,  max_rms);
    rms_post = normalize_rms(rms_post, max_rms);

    fprintf('    PRE:  M1=%.3f M2=%.3f M3=%.3f M4=%.3f\n', ...
        mean(rms_pre.exec(:,1)),  mean(rms_pre.exec(:,2)), ...
        mean(rms_pre.exec(:,3)),  mean(rms_pre.exec(:,4)));
    fprintf('    POST: M1=%.3f M2=%.3f M3=%.3f M4=%.3f\n', ...
        mean(rms_post.exec(:,1)), mean(rms_post.exec(:,2)), ...
        mean(rms_post.exec(:,3)), mean(rms_post.exec(:,4)));

    % ── iEMG ──────────────────────────────────────────────
    fprintf('  [2] iEMG...\n');
    iemg_pre  = compute_iemg(env_pre,  FS, WIN_EXEC);
    iemg_post = compute_iemg(env_post, FS, WIN_EXEC);

    fprintf('    PRE  (µV·s): M1=%.2f M2=%.2f M3=%.2f M4=%.2f\n', ...
        mean(iemg_pre(:,1)),  mean(iemg_pre(:,2)), ...
        mean(iemg_pre(:,3)),  mean(iemg_pre(:,4)));
    fprintf('    POST (µV·s): M1=%.2f M2=%.2f M3=%.2f M4=%.2f\n', ...
        mean(iemg_post(:,1)), mean(iemg_post(:,2)), ...
        mean(iemg_post(:,3)), mean(iemg_post(:,4)));

    % ── CCI ───────────────────────────────────────────────
    fprintf('  [3] CCI...\n');
    cci_pre  = compute_cci(env_pre,  FS, WIN_EXEC);
    cci_post = compute_cci(env_post, FS, WIN_EXEC);
    cci_stim = [];
    if ~isempty(env_stim)
        cci_stim = compute_cci(env_stim, FS, WIN_EXEC);
        fprintf('    STIM: mean=%.2f std=%.2f\n', mean(cci_stim), std(cci_stim));
    end
    fprintf('    PRE:  mean=%.2f std=%.2f\n', mean(cci_pre),  std(cci_pre));
    fprintf('    POST: mean=%.2f std=%.2f\n', mean(cci_post), std(cci_post));

    % ── Store ──────────────────────────────────────────────
    all_emg_results.(label).rms_pre   = rms_pre;
    all_emg_results.(label).rms_post  = rms_post;
    all_emg_results.(label).iemg_pre  = iemg_pre;
    all_emg_results.(label).iemg_post = iemg_post;
    all_emg_results.(label).cci_pre   = cci_pre;
    all_emg_results.(label).cci_post  = cci_post;
    all_emg_results.(label).cci_stim  = cci_stim;

    % ── Figures ────────────────────────────────────────────
    fig = plot_emg_envelope(env_pre, env_post, FS, TRIAL_SAMPLES, label, EMG_BIPO);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_EMG_envelope.png', label)));
    close(fig);

    fig = plot_cci(cci_pre, cci_post, cci_stim, label);
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_EMG_CCI.png', label)));
    close(fig);

    fig = figure('Visible','off','Position',[100 100 700 450]);
    x = 1:4;
    bar(x-0.2, mean(iemg_pre),  0.35, 'FaceColor',[0.2 0.4 0.8], 'DisplayName','PRE');
    hold on;
    bar(x+0.2, mean(iemg_post), 0.35, 'FaceColor',[0.8 0.2 0.2], 'DisplayName','POST');
    set(gca,'XTick',x,'XTickLabel',{'M1','M2','M3','M4'});
    ylabel('iEMG (µV·s)');
    title(sprintf('iEMG during Execution — %s', strrep(label,'_',' ')));
    legend('Location','best'); grid on;
    saveas(fig, fullfile(ROOT_RESULTS, sprintf('%s_EMG_iEMG.png', label)));
    close(fig);

    fprintf('  [OK] Figures saved\n');
end

save(fullfile(ROOT_RESULTS, 'all_emg_results.mat'), 'all_emg_results');
fprintf('\nModule 3 complete.\n');


% =========================================================
% HELPER FUNCTIONS
% =========================================================

function run_num = get_run_number(fname)
    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok); run_num = str2double(tok{1}{1});
    else; run_num = 0; end
end


function envelope = process_emg_block(file_list, folder, fs, trial_samples)
    data_all = [];
    for f = 1:length(file_list)
        EEG      = pop_loadset('filename', file_list{f}, 'filepath', folder);
        data_all = [data_all, EEG.data]; %#ok<AGROW>
    end

    bipolar = zeros(4, size(data_all,2));
    for i = 1:4
        bipolar(i,:) = data_all(2*i-1,:) - data_all(2*i,:);
    end

    nyq = fs/2;
    [b,   a  ] = butter(4, [20/nyq, 500/nyq], 'bandpass');
    [b_lp,a_lp] = butter(4, 8/nyq, 'low');

    emg_env = zeros(size(bipolar));
    for ch = 1:4
        sig = filtfilt(b, a, bipolar(ch,:));
        emg_env(ch,:) = filtfilt(b_lp, a_lp, abs(sig));
    end

    n_trials = floor(size(emg_env,2) / trial_samples);
    envelope = reshape(emg_env(:,1:n_trials*trial_samples), 4, trial_samples, n_trials);
end


function rms = compute_rms_phases(envelope, fs, win_base, win_prep, win_exec, win_rec)
    windows  = struct('baseline',win_base,'prep',win_prep,'exec',win_exec,'recovery',win_rec);
    ph_names = fieldnames(windows);
    for ph = 1:length(ph_names)
        pname = ph_names{ph};
        s = round(windows.(pname)(1)*fs)+1;
        e = round(windows.(pname)(2)*fs);
        seg = envelope(:,s:e,:);
        rms.(pname) = squeeze(sqrt(mean(seg.^2,2)))';
    end
end


function rms = normalize_rms(rms, max_val)
    for fn = fieldnames(rms)'
        rms.(fn{1}) = rms.(fn{1}) / max_val;
    end
end


function cci = compute_cci(envelope, fs, win_exec)
    s   = round(win_exec(1)*fs)+1;
    e   = round(win_exec(2)*fs);
    seg = envelope(:,s:e,:);
    n   = size(seg,3);
    cci = zeros(n,1);
    for tr = 1:n
        act   = mean(seg(:,:,tr),2);
        pair1 = 2*min(act(1),act(3))/(act(1)+act(3)+1e-10)*100;
        pair2 = 2*min(act(2),act(4))/(act(2)+act(4)+1e-10)*100;
        cci(tr) = (pair1+pair2)/2;
    end
end


function iemg = compute_iemg(envelope, fs, win_exec)
    s   = round(win_exec(1)*fs)+1;
    e   = round(win_exec(2)*fs);
    seg = envelope(:,s:e,:);
    n   = size(seg,3);
    iemg = zeros(n,4);
    for tr = 1:n
        for ch = 1:4
            iemg(tr,ch) = sum(seg(ch,:,tr))/fs;
        end
    end
end


function fig = plot_emg_envelope(env_pre, env_post, fs, trial_samples, label, ch_names)
    t   = (0:trial_samples-1)/fs;
    fig = figure('Visible','off','Position',[100 100 1200 700]);
    for ch = 1:4
        subplot(2,2,ch);
        plot(t, squeeze(mean(env_pre(ch,:,:),3)),  'b','LineWidth',1.5); hold on;
        plot(t, squeeze(mean(env_post(ch,:,:),3)), 'r','LineWidth',1.5);
        xline(3.0,'--m','Stimulus','LabelVerticalAlignment','bottom','FontSize',8,'HandleVisibility','off');
        xline(5.0,'--g','MovInit', 'LabelVerticalAlignment','bottom','FontSize',8,'HandleVisibility','off');
        xline(7.5,'--r','MovEnd',  'LabelVerticalAlignment','bottom','FontSize',8,'HandleVisibility','off');
        title(ch_names{ch},'FontSize',11);
        xlabel('Time (s)'); ylabel('Amplitude (µV)');
        legend('PRE','POST','Location','best','FontSize',8); grid on;
    end
    sgtitle(sprintf('EMG Envelope — PRE vs POST — %s', strrep(label,'_',' ')));
end


function fig = plot_cci(cci_pre, cci_post, cci_stim, label)
    fig = figure('Visible','off','Position',[100 100 700 500]);
    if ~isempty(cci_stim)
        data = {cci_pre,cci_stim,cci_post}; lbls = {'PRE','STIM','POST'};
        cols = {[0.2 0.4 0.8],[0.8 0.6 0.1],[0.8 0.2 0.2]};
    else
        data = {cci_pre,cci_post}; lbls = {'PRE','POST'};
        cols = {[0.2 0.4 0.8],[0.8 0.2 0.2]};
    end
    hold on;
    for i = 1:length(data)
        v = data{i};
        q = quantile(v,[0.25 0.5 0.75]);
        iqrv = q(3)-q(1);
        wlo  = max(min(v), q(1)-1.5*iqrv);
        whi  = min(max(v), q(3)+1.5*iqrv);
        rectangle('Position',[i-0.3,q(1),0.6,q(3)-q(1)],...
            'FaceColor',[cols{i},0.4],'EdgeColor',cols{i},'LineWidth',1.5);
        line([i-0.3,i+0.3],[q(2),q(2)],'Color',cols{i},'LineWidth',2.5);
        line([i,i],[wlo,q(1)],'Color',cols{i},'LineWidth',1.5);
        line([i,i],[q(3),whi],'Color',cols{i},'LineWidth',1.5);
        plot(i,mean(v),'o','MarkerFaceColor',cols{i},'MarkerEdgeColor','k','MarkerSize',8);
        text(i,whi+2,sprintf('%.1f±%.1f',mean(v),std(v)),...
            'HorizontalAlignment','center','FontSize',9);
    end
    set(gca,'XTick',1:length(lbls),'XTickLabel',lbls);
    ylabel('Co-Contraction Index (%)');
    title(sprintf('CCI — %s', strrep(label,'_',' '))); grid on;
end
