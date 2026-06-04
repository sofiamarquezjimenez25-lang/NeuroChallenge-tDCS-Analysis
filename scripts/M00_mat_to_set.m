% =========================================================
% MODULE 0 — Load .mat files and create EEGLAB .set files
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Reads the .mat files exported from Python and creates
% one EEGLAB .set file per run per modality (EEG and EMG).
% =========================================================

clear; clc; close all;

% ── Start EEGLAB ─────────────────────────────────────────
eeglab nogui;

% ── Paths — edit these to match your machine ─────────────
ROOT_MAT  = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\mat_files';
ROOT_SETS = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets';

% ── Constants ─────────────────────────────────────────────
FS = 2400;

EEG_CH_NAMES = {'C3','C1','Cz','C2','C4','CP1','CP2','Pz'};
EMG_CH_NAMES = {'M1L','M1R','M2L','M2R','M3L','M3R','M4L','M4R'};

% ── Session table ─────────────────────────────────────────
session_info;   % loads session_table variable

report = {};

% ── Main loop ─────────────────────────────────────────────
for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};

    mat_folder = fullfile(ROOT_MAT,  participant, session);
    set_folder = fullfile(ROOT_SETS, participant, session);

    if ~isfolder(mat_folder)
        fprintf('[SKIP] P%s S%s — mat folder not found\n', participant, session);
        continue;
    end

    if ~isfolder(set_folder)
        mkdir(set_folder);
    end

    fprintf('\n=== P%s | S%s (%s) ===\n', participant, session, sess_type);

    % Get all .mat files sorted
    mat_files = dir(fullfile(mat_folder, '*.mat'));
    names = {mat_files.name};
    [~, idx] = sort(names);
    mat_files = mat_files(idx);

    if isempty(mat_files)
        fprintf('  [WARNING] No .mat files found\n');
        continue;
    end

    for f = 1:length(mat_files)

        mat_path = fullfile(mat_folder, mat_files(f).name);

        % ── Load .mat ─────────────────────────────────────
        m = load(mat_path);

        block    = char(m.block);
        n_samp   = double(m.n_samples);

        fprintf('  [%02d] %s | block: %s | %d samples\n', ...
            f, mat_files(f).name, block, n_samp);

        % ── Build EEG dataset ─────────────────────────────
        EEG_eeg = eeg_emptyset();
        EEG_eeg.data     = double(m.eeg);       % (8, samples)
        EEG_eeg.nbchan   = 8;
        EEG_eeg.pnts     = size(m.eeg, 2);
        EEG_eeg.trials   = 1;
        EEG_eeg.srate    = FS;
        EEG_eeg.xmin     = 0;
        EEG_eeg.xmax     = (EEG_eeg.pnts - 1) / FS;
        EEG_eeg.chanlocs = make_chanlocs(EEG_CH_NAMES);
        EEG_eeg.setname  = sprintf('P%s_S%s_%s_%s_EEG', ...
            participant, session, sess_type, block);

        % ── Insert events into EEG dataset ────────────────
        EEG_eeg = insert_events(EEG_eeg, m);
        EEG_eeg = eeg_checkset(EEG_eeg);

        % ── Build EMG dataset ─────────────────────────────
        EEG_emg = eeg_emptyset();
        EEG_emg.data     = double(m.emg);       % (8, samples)
        EEG_emg.nbchan   = 8;
        EEG_emg.pnts     = size(m.emg, 2);
        EEG_emg.trials   = 1;
        EEG_emg.srate    = FS;
        EEG_emg.xmin     = 0;
        EEG_emg.xmax     = (EEG_emg.pnts - 1) / FS;
        EEG_emg.chanlocs = make_chanlocs(EMG_CH_NAMES);
        EEG_emg.setname  = sprintf('P%s_S%s_%s_%s_EMG', ...
            participant, session, sess_type, block);

        EEG_emg = insert_events(EEG_emg, m);
        EEG_emg = eeg_checkset(EEG_emg);

        % ── Save .set files ───────────────────────────────
        base_name = strrep(mat_files(f).name, '.mat', '');

        fname_eeg = sprintf('%s_EEG.set', base_name);
        fname_emg = sprintf('%s_EMG.set', base_name);

        pop_saveset(EEG_eeg, 'filename', fname_eeg, 'filepath', set_folder);
        pop_saveset(EEG_emg, 'filename', fname_emg, 'filepath', set_folder);

        report{end+1} = sprintf('P%s S%s %s %s → OK', ...
            participant, session, block, mat_files(f).name);

    end
end

% ── Final report ──────────────────────────────────────────
fprintf('\n\n========== MODULE 0 REPORT ==========\n');
for i = 1:length(report)
    fprintf('%s\n', report{i});
end
fprintf('=====================================\n');
fprintf('Done. .set files saved to: %s\n', ROOT_SETS);


% =========================================================
% HELPER FUNCTIONS
% =========================================================

function EEG = insert_events(EEG, m)
% Inserts BCI2000 events from the .mat struct into EEGLAB format

    event_fields = {
        'ev_trial_init',  'TrialInit';
        'ev_target_show', 'TargetShow';
        'ev_mov_init',    'MovInit';
        'ev_mov_end',     'MovEnd';
    };

    events = struct('type', {}, 'latency', {}, 'duration', {});

    for e = 1:size(event_fields, 1)
        field_name   = event_fields{e, 1};
        event_label  = event_fields{e, 2};

        if isfield(m, field_name)
            latencies = double(m.(field_name));  % already 1-based from Python
            for k = 1:length(latencies)
                ev.type     = event_label;
                ev.latency  = latencies(k);
                ev.duration = 1;
                events(end+1) = ev;
            end
        end
    end

    % Sort by latency
    if ~isempty(events)
        lats = [events.latency];
        [~, order] = sort(lats);
        EEG.event = events(order);
    end

end


function chanlocs = make_chanlocs(names)
% Creates minimal chanlocs struct from channel name cell array

    chanlocs = struct();
    for i = 1:length(names)
        chanlocs(i).labels = names{i};
        chanlocs(i).type   = '';
        chanlocs(i).X      = [];
        chanlocs(i).Y      = [];
        chanlocs(i).Z      = [];
    end
end