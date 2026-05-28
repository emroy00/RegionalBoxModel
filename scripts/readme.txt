5/28/26 - Eric M. Roy (emroy@mit.edu)

scripts/ contains subdirectories with all scripts required for RegionalBoxModel paper (DOI:).

figures:
- Contains four scripts:
	- function.py: Important functions required for processing data prior to plotting
	- EST_figures.ipynb: Plots all figures in main text, along with any percentages discussed therein.
	- EST_supplement.ipynb: Plots all supplemental figures.
	- EST_graphabstract.ipynb: Plots map required for background of graphical abstract figure.

preprocess:
- Contains six scripts:
	- BulkFlux.py: Generates optimal solution space for the selected year and month. Note - will only work if data is available!
        - BulkFlux_June_Supplemental: Generates optimal solution space for supplemental figures S5-S8 for selected year and month. Note - will only work if data is available in pointer directories.
        - function.py: identical to function.py in figures (soft link in original directory).
        - gc_preprocess.sh: Used to select, concatenate files from raw GC output that are then directed to data/GC_output/
        - gc_preprocess_STND.sh: Same as gc_preprocess.sh, except for STND simulation (needed to read parameters for hybrid sigma vertical coordinates used by GEOS-Chem)
        - GC_UW_FT.py: script for subsampling raw GEOS-Chem output to compile concentrations within the pbl (C_PBL), upwind (C_UW), and the free troposphere (C_FT)

