deuts is the library containing all of the relevant functions that the other programs call. It is required.

DeuteronAnalysisComponentsAllXi is the go-to program for fitting a deuterated hydrocarbon signal. Follow the commented instructions for formatting ExampleCSV1. If fitting deuterated ammonia (or other single-bond material), use DeuteronAnalysisComponentsAllXi_Single. Either way, it allows all parameters to float, outputting the found parameters, and a few plots. Use DeuteronAnalysisComponentsAllXi_Flat if the signal has already been baseline and cubic subtracted.

If this is your first time, follow these steps:
* Format ExampleCSV1.csv like so:
Line 1: Central Frequency | Frequency Span | NMR Width (number of data points)
Line 2: Baseline data
Line 3: Signal data
* Have deuts.py, ExampleCSV1,csv, and DeuteronAnalysisComponentsAllXi.py all in the same folder
* run DeuteronAnalysisComponentsAllXi.py (IDLE, Jupyter, or just run it as a file)
* check the baseline subtraction. it should look like a "batman" silhouette on a slight curve. if it doesn't, check that you put the baseline and data in ExampleCSV1 correctly
* check the output parameters - if any are "railing" (going to X.00000 or X.999999) adjust the bounds in the curve_fit until it stops railing

If the parameters are already known, use DeuteronAnalysisComponentsSimpleXi, it will run faster. You will need to edit the code to include the known parameters.

For a hole-burnt signal, use DeuteronAnalysisComponentsSimpleXi_HB (this assumes a hydrocarbon, you will need to adjust it if it's ammonia). This requires you to have fit an equilibrium signal first and gotten the parameters. The current version of this code does a wing fit on the sides, as trying to do a full cubic subtraction led to messy results. However, this means this code can only be used on fully real signals. I would recommend implementing a trial pass of the hole-burnt signal pretending it's an equilibrium signal just to get a cubic subtraction, but it will need some testing.

ExampleCSV_HB1.csv has baseline on line 2, equilibrium data on line 3, and hole-burnt data on line 4.

