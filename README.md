# Peak_Detection
This can auto detect peaks in a polymer sample in an intensity vs q plot. 


# Peak Diagnostics for 1D SAXS/WAXS Data (`peak_diag.py`)

A small utility to:
- load whitespace `.dat` files containing **q, I(q), dI(q)**
- split data into **SAXS** and **WAXS** regions by a configurable `q_split`
- detect peaks automatically (SciPy `find_peaks`) and/or overlay **manual peak positions**
- export:
  - peak tables as CSV
  - publication-ready PNG figures (via Plotly + Kaleido)
  - interactive Plotly figures (optional)

---

## Input format

By default this script expects a whitespace-delimited file with at least 3 columns:

| Column | Meaning |
|-------:|---------|
| 0 | q (Å⁻¹) |
| 1 | I(q) |
| 2 | dI(q) (uncertainty in intensity) |

Notes:
- Default behavior skips the first **2 rows** (`--skiprows 2`) to ignore header lines.
- If your file has no header, use `--skiprows 0`.

---

## Installation

### Option A: install dependencies into your current environment

```bash
pip install numpy pandas plotly scipy


## Usage + Outputs + Citation

### Usage (how to run)

Basic:
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat"

Disable saving CSV/PNG outputs:
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat" --no_save

Disable interactive display:
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat" --no_show

If your `.dat` file has no header rows:
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat" --skiprows 0

Change SAXS/WAXS split boundary (q_split in Å⁻¹):
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat" --q_split 1.2

Set wavelength (Å):
python peak_diag.py --data_dir "/path/to/data" --file "Sample_1.dat" --lambda_A 0.7863


### Outputs (where files are saved)

Outputs are saved next to your data, inside:
<data_dir>/peak_diag_output/

Example:
If you run
python peak_diag.py --data_dir "/Users/you/MyData" --file "Sample_1.dat"

Then outputs will be written to:
 /Users/you/MyData/peak_diag_output/

Typical output files:
- Sample_1_SAXS_peaks.csv
- Sample_1_WAXS_peaks.csv
- Sample_1_SAXS.png
- Sample_1_WAXS.png


### Citation

If you use this code in academic work, please cite it as:

Muhammad Faheem Hassan. "Peak Diagnostics for 1D SAXS/WAXS Data (peak_diag.py)". GitHub repository, 2026.



