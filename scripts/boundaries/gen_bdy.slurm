#!/bin/bash

#SBATCH --partition=standard
#SBATCH --qos=high
#SBATCH --mem=50G
#SBATCH --time=12:00:00                    
#SBATCH --account=class

#JVM="jre-11"
#JVM="jre-1.8.0-openjdk"
#export JAVA_HOME=/usr/lib/jvm/$JVM/
#ls /usr/lib/jvm

./iter_gen_bdy_NAARC.sh
