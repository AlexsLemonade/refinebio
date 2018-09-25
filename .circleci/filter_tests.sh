echo $(git log --format=oneline -n 1 $CIRCLE_SHA1); 
if [[ $(git log --format=oneline -n 1 $CIRCLE_SHA1) = *"noslow"* ]]; 
then 
	echo "Skipping slow tests.."; 
	./workers/run_tests.sh --exclude-tag=slow; 
else 
	echo "Running all tests.."; 
	./workers/run_tests.sh; 
fi
