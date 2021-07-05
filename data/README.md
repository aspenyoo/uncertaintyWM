# Data

This folder contains data collected from the study published in the following paper:

Yoo, A.H., Acerbi, L., & Ma, W.J. (in press). Uncertainty is Maintained and Used in Working Memory, *Journal of Vision*

The `Ellipse` folder contains data from the Ellipse condition, and the `Line` folder contains data from the Line condition. 

There are two data file types with the following naming conventions:
- condition_subjectID.csv
- subjectID_condition_simple.mat

### condition_subjectID.csv 

Column headers
- target: location of target with orientation change. 0 if a no change trial.
- response: participants' response (0: "no change", 1: "change")
- RT: reaction time, in ms
- delta: the amount of orientation change, in radians. 
- orientationI: orientation of the Ith item during the first stimulus presentation
presentation
- xI: horizontal position of Ith item on screen, in pixels, where 0 is leftmost position
- yI: vertical position of Ith item on screen, in pixels, where 0 is the topmost position
- eccI: ellipse eccentricity of the Ith item on screen

Each row is a trial.

### subjectID_condition_simple.mat
These files contain just the data used to model the data, which were published in the paper. The file contains a struct, with the following fields:
- Delta: 2000 (number of trials) x 4 (number of items) matrix, containing the orientation change of each item, in radians.
- rel: 2000 (number of trials) x 4 (number of items) matrix, containing the reliabilities (ellipse eccentricities) of each item
- resp: vector of length 2000. particiapnts response (0: "no change", 1: "change")
- pres2stimuli: string of the condition. 'Ellipse' or 'Line'
- subjid; string of the subject ID