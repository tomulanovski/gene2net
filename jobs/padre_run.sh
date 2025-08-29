#!/bin/bash
#SBATCH --job-name=padre_run
#SBATCH --output=/groups/itay_mayrose/tomulanovski/gene2net/padre_%j.out
#SBATCH --error=/groups/itay_mayrose/tomulanovski/gene2net/padre_%j.err
#SBATCH --time=100:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --partition=itaym
#SBATCH --account=itaym-users

INPUT_data="/groups/itay_mayrose/tomulanovski/gene2net/padre_test.tre"
 
java -Djava.awt.headless=true -Xmx4g -jar "/groups/itay_mayrose/tomulanovski/padre/padre-cli.jar" -i "$INPUT_data"
