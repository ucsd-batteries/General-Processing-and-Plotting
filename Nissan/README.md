# Nissan Processing and Plotting

## Setup

If you're using conda:

```bash
conda env create -f environment.yml
```

## Functionality of different scripts

- Nissan_calendar_plotting.py
  Generates SOH and IR plots for calendar aging test
- Nissan_plotting.py
  Generates SOH and IR plots for cycle aging test
- Nissan_Processing.py
  Processes raw data gathered from BMS and store them in NPx_test_summary.csv

## How to run the scripts

```bash
conda activate battery # activate conda environment

python3 Nissan_Processing.py # processes raw data into summary csv
python3 Nissan_plotting.py # plots SOH / IR for cycle aging
python3 Nissan_calendar_plotting.py # plots SOH for calendar aging
```


## Todos
- [ ] fix gantt chart
- [ ] move constants to `constants.py`
- [ ] refactor code
- [ ] automate raw data injestion
- [ ] UI for easier usage? 