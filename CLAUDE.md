# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A MATLAB code base for analyzing multi-electrode array (MEA) recordings from the
retinas of experimental animals (mice, rats, monkeys). The code loads and analyzes
data that has **already been spike-sorted and pre-analyzed by the Vision software**
(a Java application) — it works with spike times, spike-triggered averages (STAs),
electrical images (EIs), receptive-field fits, and cell classifications, not raw
voltage traces.

There is no build system. This is a library of MATLAB functions and scripts run
interactively from a MATLAB session.

## Environment setup

At the start of a MATLAB session the code base and the Vision Java library must be
on the path:

```matlab
addpath(genpath('~/Development/matlab_base/code'));   % adjust to local location
javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
```

Most `load_*` functions require the Vision `.jar` to read Vision output files, and
require the lab analysis server to be mounted (see **Data paths** below).

MATLAB is available locally, so code can be run and syntax-checked here. However,
nearly every analysis function needs live experimental data (a Vision-analyzed
`dataXXX` directory) plus the mounted server, so most functions cannot be exercised
end-to-end in this environment — reason from the source and the function headers.

## Testing

There is no comprehensive automated test suite. A copy of the `matlab_xunit`
framework is bundled at `specific_path_code/matlab_xunit_2.0.1/` and can run tests
via `runtests`. Files named `test_*.m` are mostly ad-hoc visual/verification scripts
(e.g. `code/lab/test_load_array_info.m`) rather than assertion-based unit tests.

## Core architecture: the `datarun` structure

Almost everything centers on a single struct called **`datarun`**. It is created by
`load_data` and progressively populated by `load_*` functions; analysis and plotting
functions then read from it. A typical workflow:

```matlab
datarun = load_data('2012-08-09-3/data002');  % builds the struct; loads nothing yet
datarun = load_params(datarun);                % fits, cell-type classifications
datarun = load_neurons(datarun);               % spike times -> datarun.spikes, triggers
datarun = load_sta(datarun);                   % spike-triggered averages
datarun = load_ei(datarun, 'all');             % electrical images
% then analyze/plot: get_psth, get_raster, plot_rf_fit, plot_time_course, ...
```

Key `datarun` fields you will encounter: `.names` (file paths), `.spikes` (per-cell
spike times), `.triggers`, `.cell_ids`, `.cell_types`, `.stas`, `.ei`, `.stimulus`
(parsed stimulus description, populated by `load_stim`). Cells are addressed by
Vision `cell_id`; use `get_cell_indices(datarun, cell_spec)` to convert a cell type
or id list into indices into the `datarun` arrays.

`load_data` accepts several path forms (short `piece/dataXXX`, full Vision prefix, or
an existing `datarun` to augment) — see its header for the full contract.

## Function naming conventions (house style — follow for new code)

The library uses consistent verb prefixes; match them when adding functions:

- `load_*`  — bring data from disk/Vision into the `datarun` struct
- `get_*`   — compute and return a quantity (usually without plotting)
- `plot_*`  — visualize
- `compute_*` — heavier numerical work, often called by `get_*`/`plot_*` wrappers

New functions should follow the documented header style and argument conventions in
the templates at the repo root: `code/template_wrapper_function.m`,
`code/template_low_level_function.m`. Conventions worth preserving: a `datarun`-in /
`datarun`-out signature for pipeline functions, a `cell_spec` argument parsed via
`get_cell_indices`, and optional parameters passed as a struct or name/value list
via `varargin` (many functions use `inputParser`).

## Repository layout

- `code/lab/` — **the core library.** Data loading, RF/STA analysis, plotting, and
  the `datarun` machinery. Path/machine-specific helpers live in `code/lab/paths/`.
- `code/projects/` — project-specific analyses (electrical_stim, glm, single cone,
  spike sorting, subretinal_toolbox, etc.). Larger and less general than `lab/`.
- `code/tutorials/` — teaching scripts: `matlab_tutorial.m` (general intro: loading
  data, spikes, rasters, STAs), `ds_tutorial.m` and `os_tutorial.m` (identifying
  direction- and orientation-selective RGCs). Good entry points for how the pieces
  fit together.
- `code/numerical/` — low-level geometry/math helpers.
- `code/examples/` — short example scripts for specific tasks.
- `code/external/` — **third-party, vendored. Do not modify.**
- `code/deprecated/` — **retired code kept for reference. Do not modify.**
- `utilities/` — general MATLAB helpers (workspace management, etc.), not
  retina-specific.
- `private/` — personal research code (e.g. `private/gfield/`) and other users'
  code (e.g. `private/USC/`). Editable, but it is someone's working code — change it
  only when explicitly directed, and keep changes scoped to the relevant subfolder.

## Data paths

Data-file locations resolve through `code/lab/paths/server_path.m`, which returns the
analysis server mount (e.g. `/Volumes/Analysis/`). This path is **machine-specific and
set per user** — do not hardcode a data path elsewhere or assume the value in
`server_path.m` is correct for the current machine. When `load_data` is given a
relative spec, `server_path()` is prepended. Other resolvers in `code/lab/paths/`
(`vision_path_stable.m`, `single_cone_path.m`, `stimulus_maps_path.m`, …) follow the
same per-user pattern.

## Notes for working in this code base

- The repository has accumulated duplicate/cruft files: Dropbox "conflicted copy"
  files, `.bak`/`.asv`/`.m.` files, `*_old.m`, and trailing-underscore variants
  (e.g. `plot_ei.m`, `plot_ei_.m`, `plot_ei__.m`). When editing, confirm you are on
  the canonical (usually the plain-named, most recently modified) version.
- Stimulus files: drifting-grating and similar analyses require the experiment's
  stimulus text file. Set `datarun.names.stimulus_path` before calling `load_stim`
  (see the DS/OS tutorials for the pattern).
