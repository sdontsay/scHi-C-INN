# scHi-C-INN

`scHi-C-INN` (Single-Cell Hi-C Imputation Neural Network) is a computational pipeline designed for processing and imputing single-cell Hi-C data. This project consists of three primary steps, each corresponding to a script that performs a specific task in the data processing and imputation workflow. The project aims to transform raw scHi-C data into imputed values, enhancing the utility and interpretability of single-cell genomic interactions.

## Project Workflow

### Step 1: Correlation Coefficient Calculation

- **Script Name**: `Strat_cor.r`
- **Description**: This R script calculates the correlation coefficients between individual cells in scHi-C data. It serves as the first step in the pipeline, setting the stage for subsequent analysis.
- **Usage**: 
  ```bash
  Rscript Strat_cor.r --input_dir "/path/to/input" --output_dir "/path/to/output" --stage "your_stage" --genome "hg19" --mcore 30

### Step 2: Nearest Neighbor Extraction
- ** Script Name**: Nearest_cells.py
- ** Description**: This Python script extracts the four nearest neighbors based on the correlation coefficients computed in Step 1. This step narrows down the focus to the most relevant cells for further imputation.
- ** Usage**:
  ```bash
  python Nearest_cells.py --base_dir "/path/to/base" --cell_num 16707 --genome_type "hg19" --correlation_dir "/path/to/correlation/results"

### Step 3: scHi-C Data Imputation
- ** Script Name**: scHi-C-INN.py
- ** Description**: The final Python script in the pipeline, it uses the data processed in the previous steps to impute new scHi-C data, achieving the primary goal of the project.
- ** Usage**:
  ```bash
  python scHi-C-INN.py --base_dir "/path/to/base" --cell_num 16707 --genome_type "hg19"

## Data Format
The input files for this project should be in the form of contact matrices, with values separated by a tab delimiter. Example data for demonstration purposes is available in the repository, located in the folders "Input_step1", "Input_step2", "Input_step3", and "Output_step3". Note that only chromosome 1 (chr1) data is provided as a demo.

## Example Data
The repository contains example data for each step of the pipeline. These data are intended to illustrate the expected input and output formats and to provide a basis for testing and understanding the scripts.
- ** Input_step1**: Contains raw scHi-C data for Step 1.
- ** Input_step2**: Contains data processed in Step 1, used as input for Step 2.
- ** Input_step3**: Contains data processed in Step 2, used as input for Step 3.
- ** Output_step3**: Contains the imputed scHi-C data, which is the output of Step 3.
