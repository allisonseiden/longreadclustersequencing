"""Run parent assignment algorithm on all patient IDs using phased VCFs.

Whatshap software (without indels flag)

module purge
module load python/3.5.0 py_packages/3.5
cd /hpc/users/richtf01/longreadclustersequencing/phasing
python3

"""

from PhasedData import PhasedData
import multiprocessing as mp
from utils import get_trio_df


patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
              '1-04460', '1-04537', '1-05443', '1-05673', '1-05846']
batch_i = str(2)

"""Testing

trio_df = get_trio_df()
ID = '1-05794'
patient = PhasedData(ID, trio_df, home_dir='/hpc/users/richtf01/')
whatshap_prefix = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                   'PacbioProject/IlluminaWhatshapVCFs/Batch' + batch_i +
                   '/{}/{}_chr{}_phased')

patient.illumina(whatshap_prefix)


"""


def get_pacbio_dataframes(ID):
    """Get phased de novo variants for PacBio data."""
    patient = PhasedData(ID)
    patient.pacbio()


def get_illumina_dataframes(ID):
    """Get phased de novo variants for Illumina data."""
    patient = PhasedData(ID)
    whatshap_prefix = ('/hpc/users/seidea02/www/PacbioProject/WhatshapVCFs/' +
                       '{}_illumina/{}_chr{}_phased')
    patient.illumina(whatshap_prefix)


def get_illumina_GMKF2_dataframes(ID):
    """Get phased de novo variants for Illumina data."""
    # provide trio_df if VCF IDs are not the family IDs
    trio_df = get_trio_df()
    patient = PhasedData(ID, trio_df, home_dir='/hpc/users/richtf01/')
    whatshap_prefix = ('/sc/orga/projects/chdiTrios/WGS_Combined_2017/' +
                       'PacbioProject/IlluminaWhatshapVCFs/Batch' + batch_i +
                       '/{}/{}_chr{}_phased')
    patient.illumina(whatshap_prefix)


if __name__ == '__main__':
    pool = mp.Pool(processes=5)
    # comment out the one which you're not using at the moment, only call one
    # of these at a time
    # pool.map(get_pacbio_dataframes, patientIDs)
    # pool.map(get_illumina_dataframes, patientIDs)
    pool.map(get_illumina_GMKF2_dataframes, patientIDs)
