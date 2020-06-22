#!/bin/sh

# Delete efs
ansible-playbook ~/.clap/groups/roles/namd_delete-efs.yml

# Delete security group
if [ -f "security_group.txt" ]
then
	sg_id=$(cat security_group.txt | grep group= | cut -d'=' -f2 | cut -d'"' -f1)
	ansible-playbook ~/.clap/groups/roles/namd_delete-sg.yml -e "group_id=$sg_id"
else
	echo "security_group.txt not found. Please run sh setup_aws.sh"
fi

# Delete placement group
ansible-playbook ~/.clap/groups/roles/namd_delete-pg.yml

# Clean
rm efs_ip.txt
rm security_group.txt
