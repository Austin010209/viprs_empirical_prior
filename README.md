# User Manual

## Train VIPRS

1. Use python version **under** `3.10`

2. Create `output` folder

3. Create `data` folder under the root directory

4. Download `chr_21.tar.gz` and `chr_22.tar.gz` from [Zenodo](https://zenodo.org/record/7036625#.ZEdA_-zMKrN). Unzip it and put it into `data` folder

5. Download marginal summary statistics data computed for human standing height on chromosome 21 and 22 from 337,205 White British individuals from the UK Biobank from [link](https://drive.google.com/drive/folders/1qbaGULJ3IFSW3qpOWh354EyCRoSPG05b?usp=sharing). Unzip it and put `annotations` folder and `sumstats` folder into `data` folder, but also keep the original gz file.

Then `data` folder should look like this:
```
/data
    - chr_21
    - chr_22
    - annotations
    - sumstats
```

6. Install required libraries
```
pip3 install -r requirements.txt
```

7. Run the script

- Use annotation (empirical prior VIPRS) and use chromosome 22:
```
python3 main.py --use_annot=true --chr=22
```

- Use baseline (spike and slab prior VIPRS) and use chromosome 22:
```
python3 main.py --use_annot=false --chr=22
```

## Analyse Result

1. Make sure the `output` folder contains the following files:

```
/output
    # Chromosome 21 Use annotation: False
    - fold_1_chr_21_use_annot_False.csv
    - fold_2_chr_21_use_annot_False.csv
    - fold_3_chr_21_use_annot_False.csv
    - fold_4_chr_21_use_annot_False.csv
    - fold_5_chr_21_use_annot_False.csv

    # Chromosome 21 Use annotation: True
    - fold_1_chr_21_use_annot_True.csv
    - fold_2_chr_21_use_annot_True.csv
    - fold_3_chr_21_use_annot_True.csv
    - fold_4_chr_21_use_annot_True.csv
    - fold_5_chr_21_use_annot_True.csv

    # Chromosome 22 Use annotation: False
    - fold_1_chr_22_use_annot_False.csv
    - fold_2_chr_22_use_annot_False.csv
    - fold_3_chr_22_use_annot_False.csv
    - fold_4_chr_22_use_annot_False.csv
    - fold_5_chr_22_use_annot_False.csv

    # Chromosome 22 Use annotation: True
    - fold_1_chr_22_use_annot_True.csv
    - fold_2_chr_22_use_annot_True.csv
    - fold_3_chr_22_use_annot_True.csv
    - fold_4_chr_22_use_annot_True.csv
    - fold_5_chr_22_use_annot_True.csv
```

2. Run the script

```
python3 analyse_result.py
```

3. The plots will be saved in `output` folder