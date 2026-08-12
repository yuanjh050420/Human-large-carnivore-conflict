# Human-large-carnivore-conflict
## Repository contents

This repository contains the core scripts and selected spatial data used for the *Biological Conservation* submission:

**Anthropogenic pressure and species richness predict human-large carnivore conflicts across China**

### Files

- `glmm_all_subsets.r`  
  Main GLMM workflow for fitting all candidate model subsets, ranking models by AICc, and exporting model selection, model averaging, variable importance, prediction, and diagnostic outputs.

- `glmm_top_model_r2.R`  
  Refits the top-ranked GLMM and calculates Nakagawa marginal and conditional R² values.

- `maxent.R`  
  MaxEnt workflow for aligning predictor layers, fitting the model, and generating conflict-risk prediction rasters, summary tables, and map outputs.

- `conflict layer.gdb/`  
  File geodatabase containing the spatial conflict-event layer used in the analyses.

- `spatial layers for glmm & maxent/`  
  Compressed spatial predictor layers used in GLMM and MaxEnt analyses, including carnivore richness, distance to protected areas, and poaching pressure. Other public data layers can be found in corresponding reference sites.

- `README.md`  
  Brief project overview.
