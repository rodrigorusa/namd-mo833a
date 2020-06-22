#!/bin/sh

cfg=$1
test_case=$2
stop=$3
results_folder=$4
arg=$5
value=$6

cd ${results_folder}
mkdir -p experimental_results
cd experimental_results
mkdir -p ${test_case}
cd ${test_case}
mkdir -p ${cfg}
cd ${cfg}
now="$(date +'%d-%m-%Y-%H-%M-%S')"
mkdir -p ${now}
cd ${now}

# Get cluster
clapp cluster list | grep "namd" > clusters.txt
tag=$(cut -d- -f2,3 clusters.txt)

if [ "$tag" = "$cfg" ]
then
	clapp cluster list | grep "id:" > clusters.txt
	cluster=$(cut -d: -f2 clusters.txt)

	echo "Running "${test_case}" on "${cluster}" with "${cfg}
	clapp cluster action ${cluster} namd run --extra "test_case="${test_case} "namd_args="${arg} ${value}

	echo "Fetching results from "${test_case}" on "${cluster}" with "${cfg}
	clapp cluster action ${cluster} namd fetch-results --extra "test_case="${test_case} "dest_folder="$(pwd)

	if [ "$stop" = "yes" ]
	then
		# Destroy cluster
		echo "Stopping " ${cluster}
		clapp cluster stop ${cluster}
	fi
else
	echo "cluster namd-" ${cfg} " not found"
fi

rm clusters.txt