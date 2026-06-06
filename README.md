# NeuroChallenge — tDCS Motor Study: EEG + EMG Analysis Pipeline

## Description

This repository contains the analysis scripts and processed results for a crossover tDCS motor study evaluating cortical and neuromuscular effects of anodal M1 stimulation during an 8-direction wrist motor task.

The repository does not include raw participant data, since the original recordings are not owned by the authors of this repository. Only analysis code, processed outputs, figures, and result summaries are included.

## Study Overview

The project analyzed synchronized EEG and EMG recordings from healthy adult participants during active tDCS and sham sessions. The objective was to evaluate whether anodal stimulation over the primary motor cortex (M1) modified cortical excitability and motor organization during wrist movement control.

Main signal modalities:

- EEG: central motor channels including C3, C1, Cz, C2, C4, CP1, CP2, and Pz
- EMG: four bipolar muscle pairs
- Task: 8-direction wrist motor task
- Acquisition system: BCI2000

## Repository Structure

| Folder | Description |
|--------|-------------|
| `scripts/` | MATLAB scripts used for EEG, EMG, EEG–EMG, and group-level analyses |
| `eeg results/` | EEG-related processed results and figures |
| `eeg results 2/` | Additional EEG result outputs |
| `emg results/` | EMG result figures and summaries |
| `cmc results/` | Corticomuscular coherence results |
| `erdvscci/` | ERD–CCI correlation results |
| `latency results/` | Corticomotor delay / latency results |
| `mrcp results/` | Movement-related cortical potential results |
| `group results/` | Group-level figures and summary results |
| `results/` | Additional processed outputs or general result files |

## Analysis Pipeline

| Module | Description |
|--------|-------------|
| EEG preprocessing | Filtering, common average reference, and artifact handling |
| ICA preprocessing | ICA decomposition and component inspection when applicable |
| ERD/ERS analysis | Alpha/mu and beta event-related desynchronization/synchronization across preparation, execution, and recovery phases |
| MRCP analysis | Movement-related cortical potentials aligned to EMG onset |
| EMG analysis | Integrated EMG, EMG envelope, and co-contraction index |
| CMC analysis | Beta-band corticomuscular coherence between EEG and EMG |
| ERD–CCI correlation | Trial-level relationship between cortical desynchronization and muscular co-contraction |
| Corticomotor delay | Delay between beta ERD onset and EMG onset |
| Group statistics | Descriptive and non-parametric group-level analysis |

## Main Result Categories

The processed results include:

- Resting-state spectral power analysis
- Beta ERD/ERS during the wrist motor task
- Movement-related cortical potentials
- Integrated EMG during movement execution
- Co-contraction index during movement execution
- Beta corticomuscular coherence
- Corticomotor delay
- ERD–CCI trial-level correlation
- Group-level comparison between active tDCS and sham

## Requirements

The scripts were developed for MATLAB-based analysis.

Recommended environment:

- MATLAB R2024a or compatible version
- EEGLAB
- EEGLAB plugins:
  - `firfilt`
  - `ICLabel`

Some scripts may require path adjustment depending on the local folder structure.

## Data Availability

Raw EEG and EMG recordings are not included in this repository due to data ownership and privacy considerations. The repository only contains analysis scripts and processed result outputs.

## Notes

This project was developed as part of a neuroengineering challenge focused on EEG–EMG analysis during tDCS motor modulation. Results should be interpreted as exploratory because of the small sample size, signal-quality limitations, and inter-participant variability.

## Citation

If using or adapting this repository, please cite:

Sofía Márquez Jiménez, "NeuroChallenge tDCS Motor Study: EEG–EMG Analysis Pipeline", 2026.

Code development and documentation were assisted using AI tools.
