# Uncertainty is maintained and used in working memory

This repository contains code used in the collection and analysis of the data associated with the following publication

Yoo, A.H., Acerbi, L., & Ma, W.J. (in press). Uncertainty is Maintained and Used in Working Memory, *Journal of Vision*


Here is a brief description of the folders. Each folder has a README file to further explain its contents. 

- **data/**: contains minimally altered data files (csv files) as well as data files used for model fitting (mat files) for Ellipse and Line conditions
- **experiment/**: contains functions used to collect data 
- **helper_functions/**: contains various functions used to fit, plot, and otherwise wrangle the data.
- **models/**: contains functions to fit models
- **analysis_scripts.m**: contains example code for how to collect and fit data, as well as code used to analyze model fits and generate figures. 

## `analysis_scripts.m`
This section provides some explanation of the file `analysis_scripts.m`. There are broadly three portions of the code: Collect/Fit your own data, Analyze data, and Create figures, each of which is explained a bit below. 

### Collect / Fit your own data!
This portion of the matlab script provides example code for how to run the experiment and how to fit the data. 

#### Run the experiment
Type `run_experiment` into MATLAB's command line. First, you will be prompted to give the subject's initials. You can provide any string (e.g., `S02`). Next, you will be prompted to give the experiment ID. You should input `R` for Reliability (main experimental block), `T` for Threshold block, or `P` for Practice block. Last, you will be prompted on which condition to collect data in (stimuli in second presentation). You should input `E` for Ellipse and `L` for line. Note that the low-reliability ellipse eccentricity is determined from the Threshold block, so that block is required before running the main Reliability block. 

#### Fit models to data
This section of code contains an example, with brief explanations, of how to run `find_ML_parameters`, the function which estimates the Maximum Likelihood parameter combination, given a model and data.


### Analyze data
This portion of the document contains code used to compare models, as reported in the publication. It contains code to compare individual models (using summed AICc/BIC and Bayesian Model Selection methods) and compare model families. 

### Create figures
This portion of the document contains code used to generate figures used in the publication. It contains sections that generate a plot of data alone (Figure 2), the model validation between Use and Ignore Uncertainty model (Figure 4), and the factorial model validation and AICc plots for the Line condition (Figure 5). Figure 5 can also be used to generate the factorial model validation and AICc plots for the Ellipse condition (Figure 6, in Supplementary).


