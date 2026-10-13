# Local vs Global Pupil Analysis

MATLAB code for testing whether pupil size during naturalistic free viewing is driven by *where* light falls on the retina, rather than by overall screen brightness alone.

For each video frame, the screen luminance is converted into gaze-centred (retinal) coordinates and divided into 1° × 1° patches. Each patch's luminance time series is passed through a biomechanical pupil model, and L1-regularised (LASSO) regression then finds which patches best predict the measured pupil. The regression weights form a spatial "sensitivity map" of the visual field for each participant.

> **Associated work:** [Author(s), *Title*, Year. Journal / thesis, DOI or URL — TO FILL IN]

---

## Overview of the pipeline

```
 EyeLink data + video frames
          │
          ▼
 ┌──────────────────────┐        ┌────────────────────┐
 │ 1. Preprocessing     │        │ 2. Model & analysis │
 │ p_area2mm            │        │ pupil_create_command│
 │ segmentParticipantClips ──►   │ pupil_model         │
 │ luminanceFrames      │        │ softclip            │
 │ retinalCoordinates   │        │ run_all_subjects_*  │
 │ rgb2lum              │        │ generate_final_figures
 └──────────────────────┘        └────────────────────┘
   → subjXX_retinalcoor.mat         → weight maps, R², figures, stats
```

---

## Repository contents

### 1. Preprocessing

| File | What it does |
|---|---|
| `p_area2mm.m` | Converts EyeLink pupil **area** (pixels) to pupil **diameter** (mm) using a spline fitted to artificial-eye calibration measurements (2, 4, 6, 8, 10 mm ↔ 402, 1425, 3049, 5130, 8240 px area). |
| `segmentParticipantClips.m` | For one participant, splits each movie's eye data into clips using the video frame timestamps. Blink samples and off-screen samples (gaze outside the 1280 × 1024 display) are masked, every raw sample is assigned to its nearest video frame, and gaze x, gaze y and pupil are summarised per frame (median, ignoring NaNs). Saves `subjXX_movieYY_raw.mat` (sample-level) and `subjXX_movieYY_aggData.mat` (frame-level). |
| `luminanceFrames.m` | Reads a movie file frame by frame and computes a per-pixel luminance map (cd/m²) using `rgb2lum`. Saves one `lum_frameNNNN.mat` per frame (`lum_map`, 1024 × 1280). Run once per movie. |
| `rgb2lum.m` | Converts RGB pixel values to luminance (cd/m²) with a measured gamma look-up table for the monitor. Summing the R, G and B outputs gives total pixel luminance. Option `1` (default) is the monitor type used for this dataset; look-up tables for other monitors are included. |
| `retinalCoordinates.m` | For each subject and movie, and for each frame: centres a window on the gaze position, sized at 2× the screen and padded with the measured background luminance (0.12 cd/m²), and averages luminance within 1° × 1° patches (≈41 px per degree at this setup). Also stores global (whole-screen) mean luminance, gaze position and pupil diameter. Saves `subjXX_retinalcoor.mat`. |

**Output of `retinalCoordinates.m`** (`subjXX_retinalcoor.mat`)

| Variable | Contents |
|---|---|
| `matrix{movie}` | frames × 3038 patch luminance values |
| `global_screen{movie}` | frames × 1 mean screen luminance |
| `frame_numbers{movie}` | frames × 1 video frame index |
| `gaze_px{movie}` | frames × 2 gaze position (px) |
| `pupil_mm{movie}` | frames × 1 pupil diameter (mm) |
| `metadata` | grid size (`nX` = 62, `nY` = 49), patch size, PPD, window size |

Patch index: `k = (i-1)*nY + j`, where `i` is the horizontal and `j` the vertical patch index. This means a weight vector reshapes directly into a map with `reshape(w, nY, nX)`.

Frames with missing gaze (blinks, off-screen) are stored as NaN.

### 2. Model and analysis

| File | What it does |
|---|---|
| `pupil_create_command.m` | Turns a luminance time series into a neural command for the pupil light reflex: a weighted sum of luminance and its (smoothed) rate of change, delayed by the parasympathetic latency and rectified at zero. The sympathetic command is zero (the arousal and superior-colliculus inputs are placeholders and are not used). |
| `pupil_model.m` | Biomechanical pupil model. A damped mass–spring system driven by the parasympathetic and sympathetic commands, returning a simulated pupil size over time. Based on Coutinho (2024), which builds on Fan & Yao (2011). |
| `softclip.m` | Smoothly limits the simulated pupil to a physiological range (1–7 mm) so the simulation stays stable. |
| `visual_field_pupil_control.m` | Original single-subject development script. `op = 1` regresses pupil directly on luminance; `op = 2` regresses pupil on pupil-model output. Includes diagnostic plots (e.g. optimal lag between model and data). The batch scripts below are the `op = 2` version run across participants. |
| `run_all_subjects_SPATIAL.m` | **Spatial-only model.** Intercept + 3038 patch predictors. |
| `run_all_subjects_SPATIAL-GLOBAL.m` | **Spatial + global model.** As above plus whole-screen mean luminance. |
| `run_all_subjects_GLOBAL.m` | **Global-only model.** Intercept + whole-screen mean luminance only, using the same pupil-model pipeline and regularisation for a fair comparison. |
| `generate_final_figures.m` | Loads all per-subject results and produces the figures and statistics tables (see below). |

#### How the regression works

1. **Predictors.** Each patch's luminance series is linearly interpolated over NaNs and run through `pupil_create_command` → `pupil_model`, giving one simulated pupil trace per patch.
2. **Target.** The measured pupil diameter is interpolated over NaNs, linearly detrended, and smoothed with a Gaussian kernel. The first `ig = 120` frames of each clip are discarded and all movies are concatenated.
3. **Fit.** LASSO regression (`fitrlinear`, SpaRSA solver) with λ chosen from `logspace(-5, 5, 30)` by 5-fold cross-validation (minimum MSE). The final model is refit at that λ.
4. **Output.** R² (computed on the full data at the chosen λ) and a 49 × 62 weight map per subject.

Pupil-model timing parameters used in these scripts: `dt = 0.0167 s`, parasympathetic delay `taup = 0.15 s`, sympathetic delay `taus = 0.64 s`. These were tuned per participant to account for physiological differences between participants. Parameter-tuning scripts are not included in this repository. The code here uses the fixed parameter values written in the scripts and is intended as a base for others to extend.

### Figures and statistics from `generate_final_figures.m`

| Output | Content |
|---|---|
| `Fig1` | Weight maps for six representative subjects (spanning the R² range) |
| `Fig2` | Group-average weight map; patches significant after Benjamini–Hochberg FDR (q < 0.05) |
| `Fig3` | R² distribution, per-subject R², central vs. peripheral weights, spatial-only vs. spatial + global R² |
| `Fig4` | Radial weight profile by eccentricity, with subject-level exponential-decay fits |
| `Fig5` | Reliability: split-half, spatial sign consistency, bootstrap |
| `FigS1` | R² vs. proportion of missing data |
| `FigS2` | Upper vs. lower visual field weights |
| `FigS3` | Permutation test (shuffled pupil, 100 permutations per subject) |
| `FigS4` | Spatial vs. global-only R² (from `run_global_only_model.m`) |
| `Statistics_Complete.csv`, `Statistics_Complete_v2.csv` | All summary statistics |

---

## Requirements

- MATLAB R2022a
- Statistics and Machine Learning Toolbox (`fitrlinear`, `ttest`, `normpdf`, `nansum`/`nanmean`)
- Image Processing Toolbox (`padarray`)
- Curve Fitting Toolbox (`fit`, `fittype`, Figure 4)

## Data

The eye-tracking data come from Munoz Lab: 51 participants freely viewing 10 naturalistic movie clips (EyeLink 1000, 1280 × 1024 display, 40.67 pixels per degree). Participant 18 was excluded for excessive data loss, leaving N = 50.

The scripts expect:

- `FV`: a cell array of participant structures (`ftimes`, `rawtrials`, `blinks`, `mclipnum`, `res.PPD`)
- `vidInfo`: movie information (`clipTimes`, `nFrames`)
- the movie files, for luminance extraction

## Usage

**Paths.** All scripts contain hard-coded absolute paths (`/Users/jaspreetdodd/Desktop/...`). Edit the path variables at the top of each script before running. Make sure the folder that `luminanceFrames.m` writes to is the same one `retinalCoordinates.m` reads from (`lumBasePath`).

**Run in this order:**

1. `luminanceFrames.m` — once per movie (edit `videoPath` and `outputDir` each time; the supplied version is set to movie 10).
2. `segmentParticipantClips(FV, vidInfo, participantIndex)` — once per participant.
3. `retinalCoordinates.m` — builds `subjXX_retinalcoor.mat` for all subjects (requires `FV` in the workspace).
4. Fit the models:
   - `run_all_subjects_SPATIAL.m` (spatial only)
   - `run_all_subjects_SPATIAL-GLOBAL.m` (spatial + global)
   - `run_all_subjects_GLOBAL.m` (global only; needs the no-global results to be present because it reuses their saved settings)
5. `generate_final_figures.m` — set `results_dir_noglobal`, `results_dir_global` and `figures_dir` to match where steps 4 saved results.

## Citation

If you use this code, please cite: [TO FILL IN]

Pupil model: Coutinho, J. (2024). Systems Neuroscience Models of Human Oculomotor Orienting (thesis); Fan, X. & Yao, G. (2011). Modeling transient pupillary light reflex induced by a short light flash. *IEEE Transactions on Biomedical Engineering, 58*(1), 36–42.

## License

[TO FILL IN]

## Contact

* Jaspreet Dodd - jaspreet_dodd @ sfu.ca
* Gunnar Blohm - gunnar.blohm @ queensu.ca
