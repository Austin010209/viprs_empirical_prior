# Prior type
PRIOR_TYPE_SPIKE_AND_SLAB = "PRIOR_TYPE_SPIKE_AND_SLAB"
PRIOR_TYPE_EMPIRICAL = "PRIOR_TYPE_EMPIRICAL"

# Hyper-parameters
TOTAL_FOLDS = 5
TOTAL_EPOCHS = 500

WEIGHT_INITIALIZE_MEAN = 0.0
WEIGHT_INITIALIZE_VAR = 0.001

# VAR_PER_SNP = 0.01
# in sas, var per SNP is 0.5/M, thus tau_beta = 1/(0.5/M)
# in empirical, var per SNP is is still 0.5/M, so
# 1/tau_k* = 0.5/M, tau_k* = M/0.5
# tau_beta = 1/( sum( 1/tau_k* )) approx 1/ (11/tau_k*)
# = tau_k* / 11
# tau_beta = M/0.5 / 11 = M/5.5 = 16296/5.5 = 2963
TAU_BETA_INIT_EMPIRICAL = 2963
LEARNING_RATE = 0.001

TAU_BETA_INIT = 32592
TAU_EPS_INIT_SAS = 0.01
TAU_EPS_INIT_EMPIRICAL = 1.0

PI_INIT_SAS = 0.01
PI_INIT_EMPIRICAL = 0.01

# TAU_BETA_S_INIT = 1.0
TAU_BETA_S_INIT = 11.0
MU_BETA_S_INIT = 0.0
GAMMA_S_INIT = 0.01

# Load data
DATA_DIR = "data"
OUTPUT_DIR = "output"
LD_STORE_FILES = "data/chr_{}/"
TRAIN_SUMSTATS_FILES = "data/sumstats/height/training/height_fold_{}.sumstats.gz"
TEST_SUMSTATS_FILES = "data/sumstats/height/test/height_test_fold_{}.csv.gz"
ANNOTATION_FILES = "data/annotations/baselineLD.{}.annot"
SUMSTATS_FORMAT_PLINK = "plink"
SUMSTATS_FORMAT_MAGENPY = "magenpy"
ANNOTATION_FORMAT_LDSC = "ldsc"
