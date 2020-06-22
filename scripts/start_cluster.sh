#!/bin/sh

cfg=$1

# Create temporary folder to store pub keys from nodes
mkdir -p ~/.clap/namd_keys

if [ -f "efs_ip.txt" ]
then
    efs_ip=$(cat efs_ip.txt | grep efs_ip= | cut -d'=' -f2 | cut -d'"' -f1)

	clapp cluster start namd-${cfg} --extra "efs_ip=$efs_ip"
else
	echo "efs_ip.txt not found. Please run sh setup_aws.sh"
fi

# Clean temporary directory
rm -rf ~/.clap/namd_keys
