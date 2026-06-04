% =========================================================
% MODULE 1 — EEG Preprocessing
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Steps applied to PRE and POST EEG blocks only:
%   1. Bandpass filter 1-40 Hz
%   2. Notch filter 60 Hz
%   3. Common Average Reference (CAR)
%   4. Save preprocessed .set to sets_clean folder
%
% STIM block is skipped for EEG (artifact from tDCS current)
% =========================================================

clear; clc; close all;

% ── Start EEGLAB ─────────────────────────────────────────
eeglab nogui;

% ── Paths ─────────────────────────────────────────────────
ROOT_SETS  = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets';
ROOT_CLEAN = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets_clean';

% ── Session table ─────────────────────────────────────────
session_info;

report = {};

% ── Main loop ─────────────────────────────────────────────
for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    sess_type   = session_table{s, 3};

    set_folder   = fullfile(ROOT_SETS,  participant, session);
    clean_folder = fullfile(ROOT_CLEAN, participant, session);

    if ~isfolder(set_folder)
        fprintf('[SKIP] P%s S%s — folder not found\n', participant, session);
        continue;
    end

    if ~isfolder(clean_folder)
        mkdir(clean_folder);
    end

    fprintf('\n=== P%s | S%s (%s) ===\n', participant, session, sess_type);

    % Get all EEG .set files sorted
    set_files = dir(fullfile(set_folder, '*_EEG.set'));
    names     = {set_files.name};
    [~, idx]  = sort(names);
    set_files = set_files(idx);

    if isempty(set_files)
        fprintf('  [WARNING] No EEG .set files found\n');
        continue;
    end

    for f = 1:length(set_files)

        fname = set_files(f).name;

        % ── Get run number from filename ──────────────────
        run_num = get_run_number(fname);
        block   = run_to_block(run_num);

        % ── Skip STIM block for EEG ───────────────────────
        if strcmp(block, 'stim')
            fprintf('  [SKIP] %s — STIM block excluded from EEG analysis\n', fname);
            continue;
        end

        fprintf('  [%02d] %s | R%02d → %s\n', f, fname, run_num, block);

        % ── Load .set ─────────────────────────────────────
        EEG = pop_loadset('filename', fname, 'filepath', set_folder);

        % ── Step 1: Bandpass filter 1-40 Hz ──────────────
        % Uses FIR filter via pop_eegfiltnew (firfilt plugin)
        EEG = pop_eegfiltnew(EEG, 1, 40);

        % ── Step 2: Notch filter 60 Hz ────────────────────
        EEG = pop_eegfiltnew(EEG, 59, 61, [], 1);

        % ── Step 3: Common Average Reference ─────────────
        % pop_reref with empty reference = CAR across all channels
        EEG = pop_reref(EEG, []);

        % ── Step 4: Check dataset ─────────────────────────
        EEG = eeg_checkset(EEG);

        % ── Step 5: Add block info to EEG struct ──────────
        EEG.etc.block      = block;
        EEG.etc.run_num    = run_num;
        EEG.etc.sess_type  = sess_type;
        EEG.etc.participant = participant;

        % ── Step 6: Save clean .set ───────────────────────
        out_fname = strrep(fname, '_EEG.set', ['_' block '_EEG_clean.set']);
        EEG.setname = out_fname(1:end-4);
        pop_saveset(EEG, 'filename', out_fname, 'filepath', clean_folder);

        % ── Log signal stats ──────────────────────────────
        amp_max = max(abs(EEG.data(:)));
        amp_std = std(EEG.data(:));
        fprintf('    saved | max: %.1f µV | std: %.1f µV\n', amp_max, amp_std);

        report{end+1} = sprintf('P%s S%s | %s | R%02d %s | max: %.1f | std: %.1f | OK', ...
            participant, session, sess_type, run_num, block, amp_max, amp_std);

    end
end

% ── Print report ──────────────────────────────────────────
fprintf('\n\n========== MODULE 1 REPORT ==========\n');
for i = 1:length(report)
    fprintf('%s\n', report{i});
end
fprintf('=====================================\n');
fprintf('Done. Clean files saved to: %s\n', ROOT_CLEAN);


% =========================================================
% HELPER FUNCTIONS
% =========================================================

function run_num = get_run_number(fname)
% Extracts run number from filename
% Example: 3S001R032026-03-10-11-12-35_EEG.set → 3

    tok = regexp(fname, 'R(\d+)\d{4}-', 'tokens');
    if ~isempty(tok)
        run_num = str2double(tok{1}{1});
    else
        run_num = 0;
    end
end


function block = run_to_block(run_num)
% Maps run number to experimental block label

    switch run_num
        case 1
            block = 'basal_pre';
        case 2
            block = 'stim';
        case {3, 4, 5, 6}
            block = 'pre';
        case {7, 8, 9, 10}
            block = 'post';
        case 11
            block = 'basal_post';
        otherwise
            block = 'unknown';
    end
end