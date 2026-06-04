% =========================================================
% MODULE 1b — Run ICA on all clean EEG sets
% NeuroChallenge: EEG + EMG + Kinematics + tDCS
%
% Runs ICA (runica, pca=7) on every *_EEG_clean.set
% Skips files where *_EEG_ica.set already exists
% Saves ICA weights only — no components removed
% =========================================================

clear; clc; close all;
eeglab nogui;

ROOT_CLEAN = 'C:\Users\sofia\Downloads\Neuro tDCS\Nuero Challenge\sets_clean';

session_info;

for s = 1:size(session_table, 1)

    participant = session_table{s, 1};
    session     = session_table{s, 2};
    folder      = fullfile(ROOT_CLEAN, participant, session);

    if ~isfolder(folder)
        fprintf('[SKIP] P%s S%s — folder not found\n', participant, session);
        continue;
    end

    fprintf('\n=== P%s | S%s ===\n', participant, session);

    set_files = dir(fullfile(folder, '*_EEG_clean.set'));
    [~, idx]  = sort({set_files.name});
    set_files = set_files(idx);

    for f = 1:length(set_files)
        fname     = set_files(f).name;
        ica_fname = strrep(fname, '_EEG_clean.set', '_EEG_ica.set');

        if isfile(fullfile(folder, ica_fname))
            fprintf('  [SKIP] %s — already done\n', fname);
            continue;
        end

        fprintf('  [%02d] %s\n', f, fname);

        EEG = pop_loadset('filename', fname, 'filepath', folder);
        EEG = pop_chanedit(EEG, 'lookup', which('standard_1005.elc'));
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', 7);
        EEG = pop_iclabel(EEG, 'default');
        EEG = eeg_checkset(EEG);

        EEG.setname = strrep(EEG.setname, 'clean', 'ica');
        pop_saveset(EEG, 'filename', ica_fname, 'filepath', folder);
        fprintf('    saved: %s\n', ica_fname);
    end
end

fprintf('\nModule 1b complete.\n');
