# precisionFDAsubmission
 GENeTres Submission scripts used in the Truth Challenge V2: Calling Variants from Short and Long Reads in Difficult-to-Map Regions

## Prepare Conda Environment

To reproduce the results of the submissions create a new conda environment using the requirements.txt file which contains the exact versions used during the submission.

```
conda create --name genetres --file requirements.txt
conda activate genetres
```

## Submission runs

To run the analysis simply edit the bash script containing the analysis you want to run and edit the variables at the top of the bash script. Then, just execute the bash script inside the activated conda environment. 
