# Nissan Processing and Plotting

- [Nissan Processing and Plotting](#nissan-processing-and-plotting)
  - [1. Setup](#1-setup)
  - [2. Functionality of different scripts](#2-functionality-of-different-scripts)
  - [3. Usage](#3-usage)
    - [Activate the virtual environment](#activate-the-virtual-environment)
    - [Download the raw data from google drive to the following folder](#download-the-raw-data-from-google-drive-to-the-following-folder)
    - [Running the scripts](#running-the-scripts)
  - [4. Todos](#4-todos)

## 1. Setup

If you're using conda:

```bash
conda env create -f environment.yml # creates a conda environment named "battery" & install required packages
conda activate battery # activate conda environment
which python  # check installation
```

## 2. Functionality of different scripts

- `Nissan_calendar_plotting.py` Generates SOH and IR plots for calendar aging test
- `Nissan_plotting.py` Generates SOH and IR plots for cycle aging test
- `Nissan_Processing.py` Processes raw data gathered from BMS and store them in NPx_test_summary.csv

## 3. Usage

### Activate the virtual environment

```bash
git pull # update summary data / code
conda activate battery # activate conda environment
```

### Download the raw data from google drive to the following folder

- for cycle aging: `Nissan/cycle_aging_data`
- for calendar aging: `Nissan/calendar_aging_data`

After running the scripts, the processed raw csv will be moved to corresponding subdirectory `Nissan/[calendar or cycle]_aging_data/processed`.

If any `*_plotting.py` script is executed, the plots generated will be saved to `Nissan/plots`.

### Running the scripts

In this particular sequence

If cycle aging:

```bash
cd Nissan # make sure you are in the Nissan processing sub-directory
python3 Nissan_Processing.py # processes raw data into summary csv
python3 Nissan_plotting.py # plots SOH / IR for cycle aging
```

If calendar aging:

```bash
cd Nissan # make sure you are in the Nissan processing sub-directory
python3 Nissan_Processing.py # processes raw data into summary csv
python3 Nissan_calendar_plotting.py # plots SOH for calendar aging
```

## 4. Todos

- [x] fix gantt chart
- [ ] move constants to `constants.py`
- [ ] refactor code
- [ ] automate raw data injestion with stable file structure
- [ ] automate downloading data from google drive
- [ ] UI for easier usage?
