#!/bin/bash
#zcat $1 | split -l${2} -d -a6 - ${3}/scratch/$4
echo 'pigz -dc $1 | split -l${2} -d -a6 - ${3}/scratch/$4'
pigz -dc $1 | split -l${2} -d -a6 - ${3}/scratch/$4
