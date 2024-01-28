# scHi-C-INN

`scHi-C-INN` (Single-Cell Hi-C Imputation Nearest Neighbor) is a computational pipeline designed for processing and imputing single-cell Hi-C data. This project consists of three primary steps, each corresponding to a script that performs a specific task in the data processing and imputation workflow. The project aims to transform normalized scHi-C data into imputed data, enhancing the utility and interpretability of single-cell genomic interactions.

## Prerequisites
The script requires the following to be installed:

- Python (3.x) \\
Dependencies: `NumPy`, `Pandas`

- R (4.x) \\
Dependencies: `hicrep`, `tseries`, `RJSONIO`, `stringr`, `combinat`

## Project Workflow

### Step 1: Correlation Coefficient Calculation

- **Script Name**: `Strat_cor.r`
- **Description**: This R script calculates the correlation coefficients between individual cells in scHi-C data. It serves as the first step in the pipeline, setting the stage for subsequent analysis.
- **Usage**: 
  ```bash
  Rscript Strat_cor.r --input_dir "/path/to/input" --output_dir "/path/to/output" --stage "your_stage" --genome "hg19" --mcore 30

### Step 2: Nearest Neighbor Extraction
- **Script Name**: Nearest_cells.py
- **Description**: This Python script extracts the four nearest neighbors based on the correlation coefficients computed in Step 1. This step narrows down the focus to the most relevant cells for further imputation.
- **Usage**:
  ```bash
  python Nearest_cells.py --input_dir "/path/to/input" --output_dir "/path/to/output" --cell_num 620 --label_dir "/path/to/label" --genome_type "hg19"

### Step 3: scHi-C Data Imputation
- **Script Name**: scHi-C-INN.py
- **Description**: The final Python script in the pipeline, uses the data processed in the previous steps to impute new scHi-C data, achieving the primary goal of the project.
- **Usage**:
  ```bash
  python scHi-C-INN.py --base_dir "/path/to/base" --inputs "input_dir" --outputs "scHi-C-INN" --cell_num 620 --genome_type "hg19" --correlation_dir "/path/to/input_step3"

## Data Format
The input files for this project should be in the form of contact matrices, with values separated by a tab delimiter. Example data for demonstration purposes is available in the repository, located in the folders "Input_step1", "Input_step2", "Input_step3", and "Output_step3". Note that only chromosome 1 (chr1) data is provided as a demo.

## Example Data
The repository contains example data for each step of the pipeline. These data are intended to illustrate the expected input and output formats and to provide a basis for testing and understanding the scripts.
- **Input_step1**: Contains raw scHi-C data for Step 1.
- **Input_step2**: Contains data processed in Step 1, used as input for Step 2.
- **Input_step3**: Contains data processed in Step 2, used as input for Step 3.
- **Output_step3**: Contains the imputed scHi-C data, which is the output of Step 3.
- **Demo usage**:
  ```bash
  Rscript Strat_cor.r --input_dir "/example_data/Input_step1" --output_dir "/example_data/Input_step2" --stage "BandNorm" --genome "hg19" --mcore 30
  python Nearest_cells.py --input_dir "/example_data/Input_step2" --output_dir "/example_data/Input_step3" --cell_num 620 --label_dir "/example_data/label_info.json" --genome_type "hg19"
  python scHi-C-INN.py --base_dir "/example_data" --inputs "Input_step3" --outputs "scHi-C-INN" --cell_num 620 --genome_type "hg19" --correlation_dir "/example_data/Input_step3" --num_neighbors 4
