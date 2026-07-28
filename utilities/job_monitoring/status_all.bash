#!/bin/bash

# This script is used to monitor the progress of calculations on clusters such as Shabyt, Muon, and others
# that are using SLURM for job scheduling. The script will print the current state of each job specified
# in the script (some minimal and straightforward manual adjustments are needed by the user), such as the
# current energy, basis size, number of optimization steps made, energy, excitation energy, SLURM status,
# JoBID, run time, and health status of the job. The script will print the output in a table format with
# columns separated by spaces. Depending on how things are setup, the output may be up to 140 characters
# wide.

# Full path to the directory where all individual job directories are located
JOBSDIR="/shared/home/sergiy.bubin/ecg/jobs"

# Flag to toggle between automatic and manual directory list mode (set true of false).
# In the automatic mode, the script will search for directories with the prefix specified by JOBPREFIX.
# In the manual mode, the script will use the list of directories specified in MANUAL_DIR_LIST
USE_MANUAL_DIR_LIST=true

# Manual list of directories (space-separated, only needed if USE_MANUAL_DIR_LIST=true; 
# if any item begins with "-" then a ruler in the form of a dashed line will be printed)
MANUAL_DIR_LIST=("B-_3Pe-01" "B-_1De-01" "---" "C_3Pe-01" "C_1Se-01" "C_3Se-01" "C_1Se-02" "C_3Se-02" "C_1Se-03" "---" "O_1De-01" "O_1Se-01")

# Prefix for the job directories (if not using manual list and doing auto-search, it is only needed when USE_MANUAL_DIR_LIST=false).
# For example, if "C_" is used then all directories with the names beginning with "C_" will be processed.
JOBPREFIX="C_"

# Flag to show the item number in the list of jobs (set true or false).
SHOW_LIST_ITEM_NUMBER=true

# Flag to show the eigenvalue number (set true or false). This number, as well as some other info, is extracted
# from the corresponding inout.txt file.
SHOW_EIGENVALUE_NUMBER=true

# Flag to show the number of cycles and steps made for the current basis size (set true or false).
SHOW_OPTSTP=true

# Flag to show the energy (set true or false)
SHOW_ENERGY=true

# Flag to compute the excitation energy with respect to the predefined ground state energy (set true or false).
# If the flag is set to false, the script will skip printing the column showing the excitation energy.
SHOW_EXCITATION_ENERY=false
# Predefined ground state energy (only needed if SHOW_EXCITATION_ENERY=true)
GROUND_STATE_ENERGY="-0.3784318600795520E+02"

# Flag to show the current status of job in SLURM, i.e. whether a job is currently running or pending (set true or false).
# If the flag is set to false, the script will skip printing the columns showing the SLURM status.
SHOW_JOB_SLURM_STATUS=true

# Flag to show the execution time of the job in a SLURM queue (set true or false).
SHOW_JOB_RUNTIME=true

# Flag to show the job health status, i.e. whether the calculations failed (e.g. because the energy could not be computed)
# or were terminated by SLURM/user. If the flag is set to false, the script will skip printing the columns showing the job
# health status. Otherwise, the script will print the status in the column named HEALTH and the last output file in column
# called LAST_OUTPUT.
# Note that the health status will be determined by searching through the output files in the job directories. Important: the order
# of the output files is out.txt, out1.txt, out2.txt, ..., out9.txt, out10.txt, out11.txt, out12.txt, ....
# Only the last file in this sequence will be looked at. If the end of that file contains a message about failure, the job health
# status will be set to FAILED. Otherwise, it will be set to OK.
SHOW_JOB_HEALTH_STATUS=true

# Extract short hostname and username
HOSTNAME=$(hostname -s)
USERNAME=$(whoami)

# Create a title for the screen output
TITLE="All calculations"
export TITLE+=" | "
export TITLE+=$USERNAME
export TITLE+="@"
export TITLE+=$HOSTNAME

# Construct separator lines (their length depends on how many columns are shown)
SEPLINE_DOUBLE="══════════════════"
SEPLINE_SINGLE="──────────────────"
SEPLINE_DASHED="------------------"
if [ "$SHOW_LIST_ITEM_NUMBER" = true ]; then
    export SEPLINE_DOUBLE+="════"
    export SEPLINE_SINGLE+="────"
    export SEPLINE_DASHED+="----"
fi
if [ "$SHOW_EIGENVALUE_NUMBER" = true ]; then
    export SEPLINE_DOUBLE+="═══"
    export SEPLINE_SINGLE+="───"
    export SEPLINE_DASHED+="---"
fi
if [ "$SHOW_OPTSTP" = true ]; then
    export SEPLINE_DOUBLE+="══════════"
    export SEPLINE_SINGLE+="──────────"
    export SEPLINE_DASHED+="----------"
fi
if [ "$SHOW_ENERGY" = true ]; then
    export SEPLINE_DOUBLE+="═════════════════════"
    export SEPLINE_SINGLE+="─────────────────────"
    export SEPLINE_DASHED+="---------------------"
fi
if [ "$SHOW_EXCITATION_ENERY" = true ]; then
    export SEPLINE_DOUBLE+="════════════════════"
    export SEPLINE_SINGLE+="────────────────────"
    export SEPLINE_DASHED+="--------------------"
fi
if [ "$SHOW_JOB_SLURM_STATUS" = true ]; then
    export SEPLINE_DOUBLE+="══════════════════════════════════"
    export SEPLINE_SINGLE+="──────────────────────────────────"
    export SEPLINE_DASHED+="----------------------------------"
fi
if [ "$SHOW_JOB_RUNTIME" = true ]; then
    export SEPLINE_DOUBLE+="════════════"
    export SEPLINE_SINGLE+="────────────"
    export SEPLINE_DASHED+="------------"
fi
if [ "$SHOW_JOB_HEALTH_STATUS" = true ]; then
    export SEPLINE_DOUBLE+="═══════════════════"
    export SEPLINE_SINGLE+="───────────────────"
    export SEPLINE_DASHED+="-------------------"
fi

# Build the list of directories to process
if [ "$USE_MANUAL_DIR_LIST" = true ]; then
    DIRS=()
    for DIRNAME in "${MANUAL_DIR_LIST[@]}"; do
        FULLPATH="$JOBSDIR/$DIRNAME"
        #if [ -d "$FULLPATH" ]; then
            DIRS+=("$FULLPATH")
        #fi
    done
else
    DIRS=("$JOBSDIR"/"$JOBPREFIX"*)
    # Filter only existing directories
    DIRS=("${DIRS[@]/#/}")
    DIRS=($(for d in "${DIRS[@]}"; do [ -d "$d" ] && echo "$d"; done))
fi

#Create a file with the current state of the queue
if [ "$SHOW_JOB_SLURM_STATUS" = true ]; then
    SQUEUEOUTPUT=temp_file_scr_output_squeue.txt
    squeue -o "%.6i %.6P %.12j %.19u %.2t %.11M %.8Q %.5D %.4C %.11R" -u $USERNAME > $SQUEUEOUTPUT
fi

# Loop over wanted job directories
printf "%s\n" "$SEPLINE_DOUBLE"
printf "%s\n" "$TITLE"
printf "%s\n" "$SEPLINE_SINGLE"
if [ "$SHOW_LIST_ITEM_NUMBER" = true ]; then
    printf "%-4s" "#"
fi
printf "%-12s" "STATE"
if [ "$SHOW_EIGENVALUE_NUMBER" = true ]; then
    printf " %-2s" "EV"
fi
printf " %-5s" "BASIS"
if [ "$SHOW_OPTSTP" = true ]; then
    printf " %9s" "O_CYC_STP"
fi
if [ "$SHOW_ENERGY" = true ]; then
    printf " %20s" "CURRENT_ENERGY"
fi
if [ "$SHOW_EXCITATION_ENERY" = true ]; then
    printf " %18s" "EXCITATION_ENERGY"
fi
if [ "$SHOW_JOB_SLURM_STATUS" = true ]; then
    printf " %5s %6s %8s %5s %5s" "QSTAT" "JOBID" "PARTIT" "NODES" "CPUS"
fi
if [ "$SHOW_JOB_RUNTIME" = true ]; then
    printf " %11s" "RUN_TIME"
fi
if [ "$SHOW_JOB_HEALTH_STATUS" = true ]; then
    printf " %6s %11s" "HEALTH" "LAST_OUTPUT"
fi
printf "\n"
printf "%s\n" "$SEPLINE_SINGLE"
ITEM_NUM=0
for STATEFULLDIR in "${DIRS[@]}"; do
    STATEDIR=$(basename "$STATEFULLDIR")
    # If the item begins with a minus sign, print a dashed ruler instead of
    # processing a directory (used as a visual separator in the list).
    if [[ "$STATEDIR" == -* ]]; then
        printf "%s\n" "$SEPLINE_DASHED"
        continue
    fi
    let ITEM_NUM++
    INOUTFILE=${STATEFULLDIR}/inout.txt
    CURRBASISSIZE=""
    CURRENERGY_UNFORMATTED=""
    CURRENERGY=""
    EIGVAL=""
    EXCITATION_ENERGY=""
    if [ "$SHOW_LIST_ITEM_NUMBER" = true ]; then
        printf "%-4s" "$ITEM_NUM"
    fi
    printf "%-12s" "$STATEDIR"
    # Extract the data from the job input and screen output files, if they exists
    if [[ -f $INOUTFILE ]]; then
        if [ "$SHOW_EIGENVALUE_NUMBER" = true ]; then
            EIGVAL=$(grep 'WHICH_EIGENVALUE' "$INOUTFILE" | awk '{print $2}')
            printf " %2s" "$EIGVAL"
        fi
        CURRBASISSIZE=$(grep 'BASIS_SIZE' "$INOUTFILE" | awk '{print $2}')
        printf " %5s" "$CURRBASISSIZE"
        if [ "$SHOW_OPTSTP" = true ]; then
            # The optimization cycle and step are on the line immediately following
            # the first line that contains "==========" (its exact line number varies).
            CURROPTCYCLE=$(awk '/==========/{getline; print $3; exit}' "$INOUTFILE")
            CURROPTSTEP=$(awk '/==========/{getline; print $4; exit}' "$INOUTFILE")
            CURROPT=$(echo "${CURROPTCYCLE}:${CURROPTSTEP}" | tr -d ' ')
            printf " %9s" "$CURROPT"
        fi
        CURRENERGY=""
        if [ "$SHOW_ENERGY" = true ]; then
            CURRENERGY_UNFORMATTED=$(grep 'CURRENT_ENERGY' "$INOUTFILE" | awk '{print $2}')
            CURRENERGY=$(printf "%.15f" "$CURRENERGY_UNFORMATTED")
            printf " %20s" "$CURRENERGY"
        fi
        EXCITATION_ENERGY=""
        if [ "$SHOW_EXCITATION_ENERY" = true ]; then
            # Calculate the difference using awk
            EXCITATION_ENERGY=$(awk "BEGIN { printf \"%.15f\", ($CURRENERGY_UNFORMATTED - $GROUND_STATE_ENERGY) }")
            printf " %18s" "$EXCITATION_ENERGY"
        fi
        STATUS=""
        JOBID=""
        PARTITION=""
        NODES=""
        CPUS=""
        if [ "$SHOW_JOB_SLURM_STATUS" = true ]; then
            if grep -q "$STATEDIR.*R\|R.*$STATEDIR" "$SQUEUEOUTPUT"; then
                STATUS="R"
                JOBID=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $1}')
                PARTITION=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $2}')
                NODES=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $8}')
                CPUS=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $9}')
            fi
            if grep -q "$STATEDIR.*PD\|PD.*$STATEDIR" "$SQUEUEOUTPUT"; then
                STATUS="PD"
                JOBID=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $1}')
                PARTITION=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $2}')
                NODES=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $8}')
                CPUS=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $9}')
            fi
            printf " %5s %6s %8s %5s %5s" "$STATUS" "$JOBID" "$PARTITION" "$NODES" "$CPUS"
        fi
        if [ "$SHOW_JOB_RUNTIME" = true ]; then
            RUNTIME=""
            if [ "$STATUS" = "R" ]; then
                RUNTIME=$(grep "$STATEDIR" "$SQUEUEOUTPUT" | awk '{print $6}')
            fi
            printf " %11s" "$RUNTIME"
        fi
        if [ "$SHOW_JOB_HEALTH_STATUS" = true ]; then
            HEALTH="UNKN"
            LASTOUTFILE=""
            # Search among out1.txt ... out1000.txt using the dichotomy method
            lo=1
            hi=1000
            last_found=-1
            while (( lo <= hi )); do
                mid=$(( (lo + hi) / 2 ))
                TESTFILENAME="${STATEFULLDIR}/out${mid}.txt"
                if [[ -f "$TESTFILENAME" ]]; then
                    # mid exists -> maybe there is a higher one, move right
                    last_found=$mid
                    lo=$(( mid + 1 ))
                else
                    # mid does not exist -> anything above mid also doesn't exist
                    hi=$(( mid - 1 ))
                fi
            done
            if (( last_found != -1 )); then
                HEALTH="OK"
                LASTOUTFILE="${STATEFULLDIR}/out${last_found}.txt"
            else
                # No outN.txt found; fall back to plain out.txt
                TESTFILENAME="${STATEFULLDIR}/out.txt"
                if [[ -f "$TESTFILENAME" ]]; then
                    HEALTH="OK"
                    LASTOUTFILE="$TESTFILENAME"
                fi
            fi
            LASTOUTFILE_BASENAME=""
            # Check if the last output file contains a message about failure.
            # The error takes place if the first word of the third line from the end is "Error"
            if [[ -f "$LASTOUTFILE" ]]; then
                LASTOUTFILE_BASENAME=$(basename "$LASTOUTFILE")
                if tail -1 $LASTOUTFILE | head -1 | awk '{ if ($1 == "exceeded" && $2 == "limit") exit 0; else exit 1}' ; then
                    HEALTH="FAILED"
                fi
                if tail -1 $LASTOUTFILE | head -1 | awk '{ if ($1 == "Error" && $6 == "cannot") exit 0; else exit 1}' ; then
                    HEALTH="FAILED"
                fi
            fi
            printf " %6s" "$HEALTH"
            printf " %11s" "$LASTOUTFILE_BASENAME"
        fi
    fi
    printf "\n"
done
printf "%s\n" "$SEPLINE_DOUBLE"

rm -rf $SQUEUEOUTPUT
