#!/bin/sh

set -e
set -x

if [ -f /etc/disk_added_date ] ; then
   echo "disk already added so exiting."
   exit 0
fi

sudo fdisk -u /dev/sdc <<EOF
n
p
1


t
8e
w
EOF

sudo pvcreate /dev/sdc1
sudo vgcreate VG /dev/sdc1
sudo lvcreate -n LV -l 100%VG VG
sudo mkfs.ext4 /dev/mapper/VG-LV

sudo sh -c 'date > /etc/disk_added_date'
