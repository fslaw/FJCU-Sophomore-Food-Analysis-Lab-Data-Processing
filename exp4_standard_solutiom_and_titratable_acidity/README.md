# Titration Curve Analysis

This repository contains scripts and data to analyse pH titration curves and compute titratable acidity. Files included:

- `standardization_and_titratable_acidity2_en.csv` — English CSV copy of the experimental data.
- `tc02.py` — English plotting script with propagated error bars (saved at 300 dpi).
- `tc03.py` — Same as `tc02.py` but annotates each data point with (x, y) labels with improved placement.
- `titration_curve_en.py` — original translated analysis script.

Requirements:

Install dependencies into a virtual environment (recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Run one of the scripts to generate the plots (examples):

```bash
python3 tc02.py
python3 tc03.py
```
# Titration analysis

This repository contains scripts to analyze titration data and plot titration curves.

Files:
- `tc02.py` — English-labeled analysis and plot with propagated error bars (300 dpi).
- `tc03.py` — Same as `tc02.py` but annotates each data point with (x, y) and improved label placement.
- `standardization_and_titratable_acidity2_en.csv` — English-translated CSV used by the scripts.

Requirements:
- Python 3.8+
- pandas
- numpy
- matplotlib

Usage:
1. Create a virtual environment and install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2. Run the script:

```bash
python3 tc03.py
```

The scripts save PNG files (300 dpi) in the project folder.
