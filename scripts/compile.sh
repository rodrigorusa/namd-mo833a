#!/bin/sh

if [ -f "efs_ip.txt" ]
then
	efs_ip=$(cat efs_ip.txt | grep efs_ip= | cut -d'=' -f2 | cut -d'"' -f1)

	clapp cluster start namd-compile --extra "efs_ip=$efs_ip"

	# Get cluster
	clapp cluster list | grep "namd" > clusters.txt
	cfg=$(cut -d- -f2 clusters.txt)

	if [ "$cfg" = "compile" ]
	then
		clapp cluster list | grep "id:" > clusters.txt
		cluster=$(cut -d: -f2 clusters.txt)

		# Destroy cluster
		echo "Stopping " ${cluster}
		clapp cluster stop ${cluster}
	else
		echo "cluster namd-compile not found"
	fi

	# Clean
	rm clusters.txt
else
	echo "efs_ip.txt not found. Please run sh setup_aws.sh"
fi
