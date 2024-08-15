#!/bin/bash

function run_test()
{   
    cd $folder
    echo "\n Running test $test_name:"
    Rscript $script
    if [ $? -eq 0 ]; then
        echo "Test $test_name passed."
        results+=0
    else 
        echo "Test $test_name failed."
        results+=1
    fi
    cd ..
}

results=()

# Read the input file line by line
while read line || [[ -n $line ]]; do
    # Extract the folder and file names using basename command
    folder=$(basename "$(dirname "$line")")
    script=$(basename "$line")
    # Run each test
    run_test
done < tests_to_run.txt

failed=0
n_test=${#results[@]}

for result in "${results[@]}";
    if [ result -eq 1 ]; then
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