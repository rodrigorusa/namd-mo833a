#!/bin/sh

results_folder=$1

################################ CFG-1 ###############################
# Start cluster
sh start_cluster.sh CFG-1

# Get cluster
clapp cluster list | grep "namd-" > clusters.txt
cfg=$(cut -d- -f2,3 clusters.txt)

if [ "$cfg" = "CFG-1" ]
then
	clapp cluster list | grep "id:" > clusters.txt
	cluster=$(cut -d: -f2 clusters.txt)

	# Run experiments
	for i in 1 2 3; do
		sh run.sh CFG-1 TC-1 no ${results_folder}
		sh run.sh CFG-1 TC-2 no ${results_folder}
		sh run.sh CFG-1 TC-3 no ${results_folder}
	done

	# Destroy cluster
	echo "Stopping " ${cluster}
	clapp cluster stop ${cluster}
else
	echo "cluster namd-CFG-1 not found"
fi

################################ CFG-2 ###############################
# Start cluster
sh start_cluster.sh CFG-2

# Get cluster
clapp cluster list | grep "namd-" > clusters.txt
cfg=$(cut -d- -f2,3 clusters.txt)

if [ "$cfg" = "CFG-2" ]
then
	clapp cluster list | grep "id:" > clusters.txt
	cluster=$(cut -d: -f2 clusters.txt)

	# Run experiments
	for i in 1 2 3; do
		sh run.sh CFG-2 TC-1 no ${results_folder}
		sh run.sh CFG-2 TC-2 no ${results_folder}
		sh run.sh CFG-2 TC-3 no ${results_folder}
	done

	# Destroy cluster
	echo "Stopping " ${cluster}
	clapp cluster stop ${cluster}
else
	echo "cluster namd-CFG-2 not found"
fi

################################ CFG-3 ###############################
# Start cluster
sh start_cluster.sh CFG-3

# Get cluster
clapp cluster list | grep "namd-" > clusters.txt
cfg=$(cut -d- -f2,3 clusters.txt)

if [ "$cfg" = "CFG-3" ]
then
	clapp cluster list | grep "id:" > clusters.txt
	cluster=$(cut -d: -f2 clusters.txt)

	# Run experiments
	for i in 1 2 3; do
		sh run.sh CFG-3 TC-1 no ${results_folder}
		sh run.sh CFG-3 TC-2 no ${results_folder}
		sh run.sh CFG-3 TC-3 no ${results_folder}
	done

	# Destroy cluster
	echo "Stopping " ${cluster}
	clapp cluster stop ${cluster}
else
	echo "cluster namd-CFG-3 not found"
fi

################################ CFG-4 ###############################
# Start cluster
sh start_cluster.sh CFG-4

# Get cluster
clapp cluster list | grep "namd-" > clusters.txt
cfg=$(cut -d- -f2,3 clusters.txt)

if [ "$cfg" = "CFG-4" ]
then
	clapp cluster list | grep "id:" > clusters.txt
	cluster=$(cut -d: -f2 clusters.txt)

	# Run experiments
	for i in 1 2 3; do
		sh run.sh CFG-4 TC-1 no ${results_folder}
		sh run.sh CFG-4 TC-2 no ${results_folder}
		sh run.sh CFG-4 TC-3 no ${results_folder}
	done

	# Destroy cluster
	echo "Stopping " ${cluster}
	clapp cluster stop ${cluster}
else
	echo "cluster namd-CFG-4 not found"
fi

rm clusters.txt
