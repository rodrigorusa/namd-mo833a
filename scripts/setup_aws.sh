#!/bin/sh

# Create placement group
ansible-playbook ~/.clap/groups/roles/namd_create-pg.yml

# Create security group
ansible-playbook ~/.clap/groups/roles/namd_create-sg.yml > security_group.txt
sg_id=$(cat security_group.txt | grep group= | cut -d'=' -f2 | cut -d'"' -f1)

# Create efs
ansible-playbook ~/.clap/groups/roles/namd_create-efs.yml -e "security_group=$sg_id" > efs_ip.txt
efs_ip=$(cat efs_ip.txt | grep efs_ip= | cut -d'=' -f2 | cut -d'"' -f1)
