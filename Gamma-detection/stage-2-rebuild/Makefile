data_per_script:=200
data_count:=2000
script_count:=10

all: output_merging.R

data: data_preprocessing.R
	Rscript ./data_preprocessing.R 

script: create_run_bash.sh data rank_90 rank_140 rank_150
	./create_run_bash.sh $(data_per_script)
	chmod +x ./rank_150
	chmod +x ./rank_140
	chmod +x ./rank_90

output: parallel.sh script
	./parallel.sh $(script_count)

output_merging.R: output
	Rscript ./output_merging.R

rank_90: rank_90.cpp
	g++ -o rank_90 rank_90.cpp

rank_140: rank_140.cpp
	g++ -o rank_140 rank_140.cpp

rank_150: rank_150.cpp
	g++ -o rank_150 rank_150.cpp

clean:
	rm -rf data script output order_test.h5 rank_150 rank_140 rank_90

.PHONY: all data script output rank_150 rank_140 rank_90 clean
