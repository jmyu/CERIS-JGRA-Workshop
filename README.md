
1. We suggest running the first 18 lines in ‘CERIS_Workshop.r’ to install the required R packages before the Workshop. 

2. Within each crop fold, there are three essential files for: 
•	Observations (Trait_records.txt), which stores the phenotypic values for each accession (line_code) in each environment (env_code); 
•	Genotype (Genotype.txt), which stores the marker information, if running genomics prediction to predict performance for new genotype is desired;
•	Environmental (Either Env_meta_info.txt or xxEnv_envParas_Dapx~.txt). 

3. Modify line 25 based on your setting.

3. To practice each crop, please modify line 28 & 30 correspondingly.

    Crop          Traits

    0MaizeG2F:    DTA, PH, or YLD

    1Sorghum:   	FTgdd

    2Rice:       	FTdap

    And run to the line 91, then modify line 95-96 based on the searching results.

