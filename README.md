# rms_aggregate

Aggregate the RMS results in one plot selecting the samples by conditions.

The RMSD computation for each sample must have been performed with the `rms.py` script available in the 
[RMS repository](https://github.com/njeanne/rms). The CSV results files must be gathered in directories by conditions.

The input file is a CSV (comma separated without header) which columns described respectively the condition, the path 
to the directory for the condition and the color in hexadecimal format for the condition. The CSV test data 
(`tests/inputs/inputs_conditions.csv`) looks like:

|   |   |   |
|---|---|---|
| insertions  | tests/inputs/insertions  | #fc030b |
| duplications | tests/inputs/duplications | #eb8c34 |
| WT  | tests/inputs/WT | #0303fc |

For each condition of the CSV file, the script will search for RMSD CSV files and aggregate them, then compute the 
median for each frame of the different sample. Finally the script will produce a plot of the median RMSD at each frame.

## Conda environment

A [conda](https://docs.conda.io/projects/conda/en/latest/index.html) YAML environment file is provided: 
`conda_env/rms_env.yml`. The file contains all the dependencies to run the script.

The conda environment is created using the command:
```shell script
# create the environment
conda env create -f conda_env/rms_env.yml

# activate the environment
conda activate rms
```

The usage of the script after activating the conda environment can be displayed with:
```shell script
./rms_aggregate.py -h
```

## Usage

The script can be tested with the data in the `test` directory, choosing the median method to aggregate the data:
```shell script
conda activate rms

./rms_aggregate.py --out results --md-time 10 --aggregation median \
--domain "domain ORF1" --format svg tests/inputs/inputs_conditions.csv

conda deactivate
```

The output are a line plot of the median RMSD by condition:

![median RMSD plot](doc/_static/RMSD.svg)

and an histogram of the median RMSD by condition:

![median RMSD plot](doc/_static/RMSD_histogram.svg)

The legend displays the conditions and the number of samples for each one.