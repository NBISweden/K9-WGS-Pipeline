#!/bin/sh

MOUNT_POINT=/scratch
[ ! -d $MOUNT_POINT ] && sudo mkdir $MOUNT_POINT
sudo mount /dev/mapper/VG-LV $MOUNT_POINT
sudo chown ubuntu:ubuntu $MOUNT_POINT
