#!/usr/bin/env bash

# make sure your working dictionary is the dysfunctional_detector

# create conda environment and download modules
export PATH="/home-user/yfdai/profile/download/Anaconda3/bin:$PATH"
source /home-user/yfdai/profile/download/Anaconda3/etc/profile.d/conda.sh
conda env create --file environment.yml
conda init
source ~/.bashrc 
#activating environment
conda activate dysfunctional_detector
## setting variables
conda env config vars set PATH="/home-user/yfdai/profile/download/Anaconda3/bin:$PATH"
conda env config vars set JAVA_HOME="/mnt/storage3/yfdai/download/jdk11/jdk-11.0.23"
conda env config vars set JRE_HOME="${JAVA_HOME}"
conda env config vars set CLASSPATH=".:\${JAVA_HOME}/lib:\${JRE_HOME}/lib"
conda env config vars set PATH="${JAVA_HOME}/bin:\$PATH"
conda env config vars set PYTHONPATH="/mnt/storage3/yfdai/download/script"
# re-activating environment so variable and PATH changes take effect
conda activate dysfunctional_detector
