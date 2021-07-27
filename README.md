# precisionFDAsubmission
 GENeTres Submission scripts used in the Truth Challenge V2: Calling Variants from Short and Long Reads in Difficult-to-Map Regions

## Prepare Conda Environment

To reproduce the results of our submissions just create a new conda environment using the requirements.txt file. This file contains the exact versions of tools used during the submission.

```
conda create --name genetres --file requirements.txt
conda activate genetres
```

## Submission runs

To run the analysis simply edit the bash script containing the analysis you want to run and edit the global variables at the top lines. Then, just execute the bash script inside the activated conda environment. 
