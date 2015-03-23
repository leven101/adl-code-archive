To run:
$> make
$> ./a.out -i input_file_name -t number_of_test_epochs > test_predictions_file

Example data:
 - 'wcOutput' is the Weak Classifier Output. The data is formatted as [worker/model_id, epoch_id, probability_class_0, ..., probability_class_n, hard_prediction, gold_label].
 - 'ibccOut' is the example output of 1000 test points from wcOutput above. 
 - 'perf.sh' is an example script that uses perf script to compute performance metrics. 
