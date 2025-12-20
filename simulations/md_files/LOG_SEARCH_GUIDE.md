# Log File Search Guide

Your complete reference for finding and searching through pipeline logs.

---

## Table of Contents
1. [Log File Naming Convention](#log-file-naming-convention)
2. [Finding Logs for Specific Networks](#finding-logs-for-specific-networks)
3. [Common Search Patterns](#common-search-patterns)
4. [Debugging Failed Runs](#debugging-failed-runs)
5. [Quick Reference Commands](#quick-reference-commands)

---

## Log File Naming Convention

### Location
All logs are in: `/groups/itay_mayrose/tomulanovski/gene2net/simulations/logs/`

### Naming Pattern
```
<script_name>_<config>_<job_id>_<array_task_id>.out
<script_name>_<config>_<job_id>_<array_task_id>.err
```

**Components:**
- `script_name`: prep_polyphest, run_polyphest, run_grampa, etc.
- `config`: conf_ils_low_10M, conf_ils_medium_10M, etc.
- `job_id`: SLURM job ID (e.g., 7441629)
- `array_task_id`: Task number in job array (1-105 for 21 networks × 5 replicates)

### Examples
```
run_polyphest_conf_ils_low_10M_7441629_99.out   ← Polyphest run, task 99
run_grampa_conf_ils_low_10M_p70_7435281_42.err  ← GRAMPA error log, task 42
prep_polyphest_conf_ils_low_10M_7441500_15.out  ← Polyphest prep, task 15
```

### Method-Specific Patterns

**Polyphest:**
```
prep_polyphest_<config>_<job>_<task>.out
run_polyphest_<config>_<job>_<task>.out
```

**GRAMPA:**
```
run_astral_<config>_<job>_<task>.out
run_grampa_<config>_<job>_<task>.out
```

**PADRE:**
```
prep_padre_<config>_<job>_<task>.out
run_padre_<config>_<job>_<task>.out
```

**MPSUGAR:**
```
prep_mpsugar_<config>_<job>_<task>.out
run_mpsugar_<config>_<job>_<task>.out
```

---

## Finding Logs for Specific Networks

### Network to Task ID Mapping

For **21 networks × 5 replicates = 105 tasks**, the mapping is:
```
task_id = (network_index × 5) + replicate
network_index = (task_id - 1) / 5
replicate = ((task_id - 1) % 5) + 1
```

### Network Index Order (0-indexed)
```
0:  Bendiksby_2011
1:  Koenen_2020
2:  Brysting_2007
3:  Lawrence_2016
4:  Diaz-Perez_2018   
5:  Wisecaver_2023
6:  Ding_2023
7:  Liang_2019
8:  Popp_2005
9:  Wu_2015
10: Liu_2023
11: Ren_2024
12: Marcussen_2011
13: Marcussen_2012
14: Sessa_2012b
15: Zhao_2021
16: Hori_2014
17: Marcussen_2015
18: Shahrestani_2015
19: Morales-Briones_2021
20: Soza_2014
```

### Calculate Task IDs for a Network

**Diaz-Perez_2018 = index 4:**
```
Replicate 1: task_id = 4×5 + 1 = 21
Replicate 2: task_id = 4×5 + 2 = 22
Replicate 3: task_id = 4×5 + 3 = 23
Replicate 4: task_id = 4×5 + 4 = 24
Replicate 5: task_id = 4×5 + 5 = 25
```

**So Diaz-Perez_2018 is in tasks 21-25!**

---

## Finding Diaz-Perez_2018 Logs

### All Polyphest Logs for Diaz-Perez_2018

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs

# Find all Polyphest logs for tasks 21-25
ls -lh run_polyphest_conf_ils_low_10M_*_{21..25}.out

# Example output:
# run_polyphest_conf_ils_low_10M_7441629_21.out
# run_polyphest_conf_ils_low_10M_7441629_22.out
# run_polyphest_conf_ils_low_10M_7441629_23.out
# run_polyphest_conf_ils_low_10M_7441629_24.out
# run_polyphest_conf_ils_low_10M_7441629_25.out
```

### Check What Network Each Log Is For

```bash
# Verify which network is in each log
for task in {21..25}; do
  echo "=== Task $task ==="
  grep "Network:" run_polyphest_conf_ils_low_10M_*_${task}.out | head -1
  echo ""
done

# Should show:
# === Task 21 ===
# Network: Diaz-Perez_2018 (index: 4)
# ...
```

### Check Replicate Numbers

```bash
# Check replicate for each task
for task in {21..25}; do
  echo "Task $task:"
  grep "Replicate:" run_polyphest_conf_ils_low_10M_*_${task}.out | head -1
done
```

### View Specific Log

```bash
# Look at Diaz-Perez_2018, replicate 1 (task 21)
less run_polyphest_conf_ils_low_10M_*_21.out

# Or view the whole log
cat run_polyphest_conf_ils_low_10M_*_21.out

# Just the end (to see if it completed)
tail -50 run_polyphest_conf_ils_low_10M_*_21.out
```

---

## Common Search Patterns

### 1. Check if Jobs Completed Successfully

```bash
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs

# Check for success messages
grep "COMPLETED SUCCESSFULLY" run_polyphest_conf_ils_low_10M_*_{21..25}.out

# Check for failure messages
grep "FAILED" run_polyphest_conf_ils_low_10M_*_{21..25}.out
```

### 2. Find Errors

```bash
# Look for ERROR in output logs
grep -i "error" run_polyphest_conf_ils_low_10M_*_{21..25}.out

# Look in error logs (.err files)
cat run_polyphest_conf_ils_low_10M_*_{21..25}.err

# Find any non-empty error logs
for task in {21..25}; do
  err_file=$(ls run_polyphest_conf_ils_low_10M_*_${task}.err 2>/dev/null)
  if [ -s "$err_file" ]; then
    echo "=== $err_file has errors ==="
    cat "$err_file"
    echo ""
  fi
done
```

### 3. Check Runtime/Duration

```bash
# See how long each job took
for task in {21..25}; do
  echo "Task $task:"
  grep "Duration:" run_polyphest_conf_ils_low_10M_*_${task}.out
done
```

### 4. Find Specific Method Runs

```bash
# All Polyphest runs for Diaz-Perez
ls run_polyphest_*_{21..25}.out

# All GRAMPA runs for Diaz-Perez
ls run_grampa_*_{21..25}.out

# All PADRE runs for Diaz-Perez
ls run_padre_*_{21..25}.out
```

### 5. Check All Methods for One Network

```bash
# Create a summary for Diaz-Perez_2018 (tasks 21-25)
echo "Diaz-Perez_2018 Pipeline Status"
echo "================================"
echo ""

for method in polyphest grampa padre mpsugar; do
  echo "$method:"
  for task in {21..25}; do
    if ls run_${method}_*_${task}.out 2>/dev/null | grep -q .; then
      status=$(grep -l "COMPLETED SUCCESSFULLY" run_${method}_*_${task}.out 2>/dev/null)
      if [ -n "$status" ]; then
        echo "  Task $task: SUCCESS"
      else
        echo "  Task $task: FAILED or RUNNING"
      fi
    else
      echo "  Task $task: NO LOG FOUND"
    fi
  done
  echo ""
done
```

---

## Debugging Failed Runs

### Step-by-Step for Diaz-Perez_2018 Polyphest

#### 1. Check if logs exist
```bash
ls -lh run_polyphest_conf_ils_low_10M_*_{21..25}.out
```

**If NO logs:** Job never ran or wrong config name

#### 2. Check job completion status
```bash
# For each replicate (tasks 21-25)
for task in {21..25}; do
  echo "=== Replicate $((task - 20)) (task $task) ==="

  # Find the log file
  log_file=$(ls run_polyphest_conf_ils_low_10M_*_${task}.out 2>/dev/null | head -1)

  if [ -z "$log_file" ]; then
    echo "NO LOG FILE FOUND"
  else
    echo "Log: $log_file"

    # Check completion
    if grep -q "COMPLETED SUCCESSFULLY" "$log_file"; then
      echo "Status: SUCCESS"
    elif grep -q "FAILED" "$log_file"; then
      echo "Status: FAILED"
    else
      echo "Status: UNKNOWN (check manually)"
    fi

    # Show last 10 lines
    echo "Last lines:"
    tail -10 "$log_file"
  fi
  echo ""
done
```

#### 3. Look for specific errors
```bash
# Common Polyphest errors
for task in {21..25}; do
  log=$(ls run_polyphest_conf_ils_low_10M_*_${task}.out 2>/dev/null | head -1)
  if [ -n "$log" ]; then
    echo "=== Task $task ==="

    # Check for file not found errors
    grep -i "not found" "$log" || echo "No 'not found' errors"

    # Check for Python errors
    grep -i "traceback" "$log" || echo "No Python tracebacks"

    # Check for timeout
    grep -i "timeout\|killed" "$log" || echo "No timeout/killed"

    echo ""
  fi
done
```

#### 4. Check prep step completed
```bash
# Polyphest needs prep before run
# Check if prep completed for Diaz-Perez

echo "Checking Polyphest PREP for Diaz-Perez_2018..."
for task in {21..25}; do
  prep_log=$(ls prep_polyphest_conf_ils_low_10M_*_${task}.out 2>/dev/null | head -1)
  if [ -n "$prep_log" ]; then
    if grep -q "completed successfully\|COMPLETED" "$prep_log"; then
      echo "Task $task prep: OK"
    else
      echo "Task $task prep: FAILED - check $prep_log"
    fi
  else
    echo "Task $task prep: NO LOG"
  fi
done
```

#### 5. Check input files were created
```bash
# Verify prep created the input files
network="Diaz-Perez_2018"
config="conf_ils_low_10M"

for rep in {1..5}; do
  input_dir="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations/${network}/processed/${config}/polyphest_input/replicate_${rep}"

  echo "Replicate $rep:"
  if [ -d "$input_dir" ]; then
    ls -lh "$input_dir"
  else
    echo "  Directory not found: $input_dir"
  fi
  echo ""
done
```

---

## Quick Reference Commands

### Universal Log Search Template

```bash
# Find logs for ANY network
cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs

# 1. Calculate task IDs
# network_index = <find in list above>
# tasks = network_index × 5 + 1 through network_index × 5 + 5

# 2. List logs
ls <script_name>_<config>_*_{task1..task5}.out

# 3. Check status
grep "COMPLETED SUCCESSFULLY\|FAILED" <script_name>_<config>_*_{task1..task5}.out

# 4. View errors
cat <script_name>_<config>_*_{task1..task5}.err
```

### Method-Specific Quick Checks

```bash
# Polyphest for network index N
tasks=$((N*5+1))
taskend=$((N*5+5))
ls run_polyphest_conf_ils_low_10M_*_{$tasks..$taskend}.out

# GRAMPA for network index N
ls run_grampa_conf_ils_low_10M_*_{$tasks..$taskend}.out

# PADRE for network index N
ls run_padre_conf_ils_low_10M_*_{$tasks..$taskend}.out
```

### Search Within Logs

```bash
# Find all logs mentioning a specific error
grep -l "division by zero" run_polyphest_*.out

# Search for network name across all logs
grep -l "Diaz-Perez_2018" run_*.out | sort

# Find logs that didn't complete
grep -L "COMPLETED SUCCESSFULLY" run_polyphest_conf_ils_low_10M_*.out

# Count successes
grep -l "COMPLETED SUCCESSFULLY" run_polyphest_conf_ils_low_10M_*.out | wc -l
```

---

## Network Index Quick Reference

For quick lookups:
```bash
# Create a helper function
get_tasks() {
  local networks=(
    "Bendiksby_2011"      # 0:  tasks 1-5
    "Koenen_2020"         # 1:  tasks 6-10
    "Brysting_2007"       # 2:  tasks 11-15
    "Lawrence_2016"       # 3:  tasks 16-20
    "Diaz-Perez_2018"     # 4:  tasks 21-25
    "Wisecaver_2023"      # 5:  tasks 26-30
    "Ding_2023"           # 6:  tasks 31-35
    "Liang_2019"          # 7:  tasks 36-40
    "Popp_2005"           # 8:  tasks 41-45
    "Wu_2015"             # 9:  tasks 46-50
    "Liu_2023"            # 10: tasks 51-55
    "Ren_2024"            # 11: tasks 56-60
    "Marcussen_2011"      # 12: tasks 61-65
    "Marcussen_2012"      # 13: tasks 66-70
    "Sessa_2012b"         # 14: tasks 71-75
    "Zhao_2021"           # 15: tasks 76-80
    "Hori_2014"           # 16: tasks 81-85
    "Marcussen_2015"      # 17: tasks 86-90
    "Shahrestani_2015"    # 18: tasks 91-95
    "Morales-Briones_2021" # 19: tasks 96-100
    "Soza_2014"           # 20: tasks 101-105
  )

  local network="$1"
  for i in "${!networks[@]}"; do
    if [[ "${networks[$i]}" == "$network" ]]; then
      local start=$((i * 5 + 1))
      local end=$((i * 5 + 5))
      echo "Network: $network (index $i)"
      echo "Tasks: $start-$end"
      return
    fi
  done
  echo "Network not found: $network"
}

# Usage:
# get_tasks "Diaz-Perez_2018"
```

---

## Diaz-Perez_2018 Specific Investigation Script

Save this as `check_diaz_perez.sh`:

```bash
#!/bin/bash

cd /groups/itay_mayrose/tomulanovski/gene2net/simulations/logs

network="Diaz-Perez_2018"
config="conf_ils_low_10M"
tasks=(21 22 23 24 25)

echo "========================================="
echo "Investigating $network ($config)"
echo "Tasks: 21-25 (5 replicates)"
echo "========================================="
echo ""

# Check each method
for method in polyphest grampa padre mpsugar; do
  echo "--- $method ---"

  found=false
  for task in "${tasks[@]}"; do
    rep=$((task - 20))

    # Find log files (handle different naming patterns)
    if [ "$method" = "polyphest" ]; then
      prep_log=$(ls prep_${method}_${config}_*_${task}.out 2>/dev/null | head -1)
      run_log=$(ls run_${method}_${config}_*_${task}.out 2>/dev/null | head -1)
    else
      run_log=$(ls run_${method}_${config}_*_${task}.out 2>/dev/null | head -1)
    fi

    if [ -n "$run_log" ]; then
      found=true

      # Check status
      if grep -q "COMPLETED SUCCESSFULLY\|completed successfully" "$run_log" 2>/dev/null; then
        status="SUCCESS"
      elif grep -q "FAILED" "$run_log" 2>/dev/null; then
        status="FAILED"
      else
        status="UNKNOWN"
      fi

      echo "  Rep $rep (task $task): $status"

      # If failed, show why
      if [ "$status" != "SUCCESS" ]; then
        echo "    Last 5 lines:"
        tail -5 "$run_log" | sed 's/^/      /'
      fi
    fi
  done

  if [ "$found" = false ]; then
    echo "  No logs found (never ran?)"
  fi

  echo ""
done

echo "========================================="
echo "Check output files:"
echo "========================================="

base_dir="/groups/itay_mayrose/tomulanovski/gene2net/simulations/simulations"

for method in grampa polyphest_p50 polyphest_p70 polyphest_p90 padre mpsugar; do
  echo "$method:"
  for rep in {1..5}; do
    if [ "$method" = "grampa" ]; then
      result_file="${base_dir}/${network}/results/${config}/${method}/replicate_${rep}/grampa_result.tre"
    elif [[ "$method" == polyphest* ]]; then
      result_file="${base_dir}/${network}/results/${config}/${method}/replicate_${rep}/polyphest_result.tre"
    elif [ "$method" = "padre" ]; then
      result_file="${base_dir}/${network}/results/${config}/${method}/replicate_${rep}/padre_result.tre"
    elif [ "$method" = "mpsugar" ]; then
      result_file="${base_dir}/${network}/results/${config}/${method}/replicate_${rep}/mpsugar_result.tre"
    fi

    if [ -f "$result_file" ]; then
      size=$(stat -c%s "$result_file" 2>/dev/null || stat -f%z "$result_file" 2>/dev/null)
      echo "  Rep $rep: EXISTS ($size bytes)"
    else
      echo "  Rep $rep: NOT FOUND"
    fi
  done
  echo ""
done
```

Run it:
```bash
chmod +x check_diaz_perez.sh
./check_diaz_perez.sh
```

This will give you a complete picture of what happened with Diaz-Perez_2018!
