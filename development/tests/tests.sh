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

# Read the input file line by line
while read line || [[ -n $line ]]; do
    # Extract the folder and file names using basename command
    folder=$(basename "$(dirname "$line")")
    script=$(basename "$line")
    test_name=$(basename "$line" .R)
    # Run each test
    run_test
done < tests_to_run.txt

failed=0

for result in "${results[@]}"; do
    if [ ! $result -eq 0 ]; then
        ((failed++))
    fi
done


if [ $failed -gt 0 ]; then
    echo "$failed out of $n_test tests failed."
    return 1
else 
    echo "All tests passed!"
    return 0
fi