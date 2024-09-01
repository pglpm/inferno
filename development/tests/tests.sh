#!/bin/bash

function run_test()
{   
    cd $folder
    ((n_test++))
    echo "------------------"
    echo "Running test # $n_test: $test_name"
    Rscript $script 2>&1 | sed "s/^/[R: $test_name] /"
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "Test $test_name passed."
        results+=(0)
    else 
        echo "Test $test_name failed."
        results+=(1)
    fi
    echo "------------------"
    cd ..
}

results=()
n_test=0
failed=0
filename=$1

# Read input argument
if [[ -f $filename ]] ; then
    echo "Running tests in $filename."
    else
        echo "Please input a valid text file containing a list of tests to run."
    fi

# Read the input file line by line
while read line || [[ -n $line ]]; do
    # Silently skip commented or empty lines 
    if [[ ${line:0:1} != "#" ]] && [[ ${line:0:1} != "" ]] ; then
    
        # Extract the folder and file names using basename command
        echo "Parsing line $line."
        folder=$(basename "$(dirname "$line")")
        script=$(basename "$line")
        test_name=$(basename "$line" .R)

        # Check if folder and file exists
        if [[ -d $folder ]] ; then
            if [[ -f $line ]] && [ "$line" != "" ]; then
                run_test
            else
                echo "$line is not a file, skipping line.";
            fi
        else
            echo "$folder is not a directory, skipping line.";
        fi
    fi
    
done < "$filename"

for result in "${results[@]}"; do
    if [ ! $result -eq 0 ]; then
        ((failed++))
    fi
done

if [ $n_test -gt 0 ]; then
    if [ $failed -gt 0 ]; then
        echo "$failed out of $n_test tests failed."
        return 1
    else 
        echo "All tests passed!"
        return 0
    fi
else
    echo "No tests were run."
    return 0
fi