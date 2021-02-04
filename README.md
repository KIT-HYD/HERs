# HER addaptation for categorical data and probability maps: 
HER method,including sequential simulation (HERs).

## License agreement

The HER method comes with ABSOLUTELY NO WARRANTY. You are welcome to modify and redistribute it within the license agreement. The HER method is published under the CreativeCommons "CC-BY-4.0" license together with a ready-to-use sample data set. To view a full version of the license agreement please visit [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

## Requisites

* MATLAB (tested on 2018b).

## Usage

See HER.m

## File structure

* HER script ... .m
* functions/ ... .m
* datasets/ ... .mat

### HER

The main (sHER_E_0X_Y_Jura_struct.m) script is divided in the following sections:

__1. Load dataset__
	Loads the dataset.
	
__2. Define infogram and HER3 properties__
	Definition of the infogram properties, aggregation method, z threshold (optional).
	
__3. HER1: Spatial characterization__
	Extracts spatial correlation patterns. ```f_her_infogram.m```

__4. HER2: Weight optimization__
	Optimizes weights for the aggregation method based on entropy minimization. ```f_her_weight.m```

__5. HER3: z PMF prediction__
	Applies spatial characterization and optimal weights for PMF prediction. ```f_her_predict.m```
	
__6. HERs4: Sequential Simulation__
Applies spatial characterization and optimal weights for PMF prediction. ```f_her_predict.m```

The sHER_E_0304_Extract_pmf_statistics_performance_and_plot.m script extract PMF statistics, calculate performance, and plot maps. It is composed by:

__1. Extract PMF statistics__
	Obtains mean, median, mode and probability of a z threshold (optional) of the predicted z PMFs and plots the results.

__2. Calculate performance metrics__
	Calculates Root Mean Square Error (RMSE), Mean Error (ME), Mean Absolute Error (MAE), Nash-Sutcliffe model efficiency and scoring rule (DKL) of the validation set.  
	
__3. Plot maps


### Functions

The functions are detailed in their own source code body. Examples of how to use them are available in the `sHER.m` script. 

```
f_plot_probabilitymap.m
```

```f_plot``` functions were specifically built for the dataset of the study.

### Dataset of the study

Each dataset file contains:

* __idx_rand_full:__ index of the randomly shuffled data (same for all files)
* __sample_size:__ all calibration sizes available of the dataset (same for all files)
* __data:__ matrix with z values of the full generated dataset
* __txt:__ dataset type (SR0, SR1, LR0, LR1)
* __idx_cal:__ index of the calibration set
* __idx_val:__ index of the validation set
* __idx_test:__ index of the test set
* __x:__ matrix with x coordinates of the full dataset
* __x_cal:__ vector with x coordinates of the calibration set (x_cal=x(idx_cal))
* __x_val:__ vector with x coordinates of the validation set (x_val=x(idx_val))
* __x_test:__ vector with x coordinates of the test set (x_test=x(idx_test))
* __y:__ matrix with y coordinates of the full dataset
* __y_cal:__ vector with y coordinates of the calibration set (y_cal=y(idx_cal))
* __y_val:__ vector with y coordinates of the validation set (y_val=y(idx_val))
* __y_test:__ vector with y coordinates of the test set (y_test=y(idx_test))
* __z:__ matrix with z values of the full generated dataset (z=data)
* __z_cal:__ vector with z values of the calibration dataset (z_cal=z(idx_cal))
* __z_val:__ vector with z values of the validation dataset (z_val=z(idx_val))
* __z_test:__ vector with z values of the test dataset (z_test=z(idx_test))
* __dim_cal:__ size of the calibration set (dim_cal=length(idx_cal))
* __dim_val:__ size of the validation set (dim_val=length(idx_val))
* __dim_test:__ size of the test set (dim_test=length(idx_test))


## Contact

Stephanie Thiesen | stephanie.thiesen@kit.edu
Uwe Ehret | uwe.ehret@kit.edu

