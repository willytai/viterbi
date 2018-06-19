====== File description ======

+- dsp_hw1/
   +-  c_cpp/
   |     +-                    some hmm program
   +-  modellist.txt           model list to train
   +-  model_init.txt          initial model for training
   +-  seq_model_01~05.txt     training data observation
   +-  testing_data1.txt       testing data  observation
   +-  testing_answer.txt      answer for "testing_data1.txt"
   +-  testing_data2.txt       testing data without answer

====== Program Execute ======

c/cpp:
 make
 ./train $iter model_init.txt seq_model_01.txt model_01.txt
 ./train $iter model_init.txt seq_model_02.txt model_02.txt
 ./train $iter model_init.txt seq_model_03.txt model_03.txt
 ./train $iter model_init.txt seq_model_04.txt model_04.txt
 ./train $iter model_init.txt seq_model_05.txt model_05.txt
 ./test modellist.txt testing_data1.txt result1.txt
 ./test modellist.txt testing_data2.txt result2.txt

$iter is positive integer, which is iteration of Baum-Welch algorithm.
