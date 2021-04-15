# Uncertainty is maintained and used in working memory

This repository contains code used in the collection and analysis of the data associated with the following publication

Yoo, A.H., Acerbi, L., & Ma, W.J. (in press). Uncertainty is Maintained and Used in Working Memory, *Journal of Vision*


Here is a brief description of the folders. Each folder has a README file to further explain its contents. 

- **data**: contains minimally altered data files (csv files) as well as data files used for model fitting (mat files) for Ellipse and Line conditions
- **experiment**: contains code used to collect data 
- **helper_functions**: various functions necessary to fit, plot, or otherwise wrangle the data.
- **models**: code to fit and compare models
- **analysis_scripts.m**: code used to collect, analyze, and plot data. 

## `analysis_scripts.mat`

### Run the experiment
Type `run_experiment` into MATLAB's command line. First, you will be prompted to give the subject's initials. You can provide any string (e.g., `S02`). Next, you will be prompted to give the experiment ID. You should input `R` for Reliability (main experimental block), `T` for Threshold block, or `P` for Practice block. Last, you will be prompted on which condition to collect data in (stimuli in second presentation). You should input `E` for Ellipse and `L` for line. Note that the low-reliability ellipse eccentricity is determined from the Threshold block, so that block is required before running the main Reliability block. 

### Fit models to data
To fit the models from this publication, you can run `find_ML_parameters` with the appropriate inputs. An example of how to run this data for a particular participant can be found in `analysis_scripts.m` in the section titled "fit data"

### Analyze data
This portion of the document contains code to compare individual models (using summed AICc/BIC and Bayesian Model Selection methods) and compare model families. 

### Create figures
This portion of the document contains code used to generate Figure 2 (data plot), 4 (model comparison between Use and Ignore Uncertainty model), and 5 (factorial model comparison) in the manuscript. 


