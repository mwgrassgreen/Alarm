# Alarm

Alarm is a R package of a wearable alarming system based on CuSum statistics. 


## Dependence
* [R](https://www.r-project.org/) (version >= 3.3.0)

* R package: "xts"

* R function file: "online_cusum_alarm_fn.R"

## Usage
To set the working directory as "/Alarm/R"

`source("online_cusum_alarm_fn.R")`

`online.alarming.fn(peo.id.1, dir.hr, dir.step, track.par=12, gap.thres = 14)`

--peo.id.1 individual ID

--dir.hr directory of raw heart rate data

--dir.step directory of raw step data

--track.par tracking parameter (default: 12 hours); tracking the trend of CuSum statistics once the initial alarm is called

--gap.thres missing data gap (default: 14 days); if the missing data gap is greater than 14 days, to restart the alarm system

The evaluation is processed through chunks of data (2 * 28 days). 


## Output
* reformatted HR data, step data, and smoothed RHR data are saved in subfolder /output/clean

* online result (based on CuSum) and offline result (based on RHR-diff) for each chunk are saved in subfolder /output/result

* online summary table, offline summary table, and daily evaluation table for all the data are saved in subfolder /output/table

* notice that the CuSum alarm cannot run due to lack of enough data for some some is saved in subfolder /output/note

* online alarm figure combined with offline detection figure for each chunk is saved in subfolder /output/figure






