# Contaminant Analysis for heavy ion accelerators (Linux version)
This code calculates possible ion contaminants based on the A/q difference. 
The charge state distribution after stripper is calculated with ETACHA4.
Beam tracking and beam loss analysis are then implemented with JuTrack (pyJuTrack).

## Setup (Required Before Running)

### 1. Configure ETACHA paths
Open `etacha_batch_run.py` and update the following variables:

- `ETACHA_EXE`: Full path to your `ETACHA4.exe` in Windows. In Linux, use `etacha4` directly.
- `RESULTS_DIR`: Directory where ETACHA stores output files  
  (typically `~/Documents/LISEcute/results`)

ETACHA4 is included with LISE and can be downloaded at:  
https://lise.frib.msu.edu/download/download.html

### 2. Provide a template input file
Copy the file `template.etacha` into the working directory.  
This file is used as the base configuration for all batch ETACHA runs.

### 3. [Optional] Customize the contaminant species range
The list of species to scan is defined as `DEFAULT_SPECIES` in  
`contaminant_calculator.py`. Edit the definition of `DEFAULT_SPECIES` if you need to restrict or expand
the range of contaminant species to search.

### 4. Install pyJuTrack for beam tracking
Beam tracking features require the Python interface to JuTrack  
(**Python 3.13 or later is required**).

Install pyJuTrack in editable mode in you Python environment:

```bash
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl.git
cd JuTrack.jl/python_integration
pip install -e .
```

## Usage

See `run_analysis.py`. The first run of JuTrack will install all dependencies that takes a few minutes.
Don't forget to clean up the result folder `RESULTS_DIR` after running the code.

## Visualize Charge Distributions before and after stripper

```python
from visualize_charge_distribution import plot_charge_distribution
contaminants_file = "my_analysis_list.csv"
charge_dist_file = "my_analysis_charge_distributions.csv"
plot_charge_distribution(contaminants_file, charge_dist_file, "charge_states.png")
```
