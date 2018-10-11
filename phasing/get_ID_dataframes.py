"""Run parent assignment algorithm on all patient IDs using phased VCFs.

Whatshap software (without indels flag)
"""

from PhasedData import PhasedData
import multiprocessing as mp

patientIDs = ['1-00801', '1-01019', '1-03897', '1-04190', '1-04389',
              '1-04460', '1-04537', '1-05443', '1-05673', '1-05846']


def get_pacbio_dataframes(ID):
    """Get phased de novo variants for PacBio data."""
    patient = PhasedData(ID)
    patient.pacbio()


def get_illumina_dataframes(ID):
    """Get phased de novo variants for Illumina data."""
    patient = PhasedData(ID)
    patient.illumina()


if __name__ == '__main__':
    pool = mp.Pool(processes=5)
    # comment out the one which you're not using at the moment, only call one
    # of these at a time
    pool.map(get_pacbio_dataframes, patientIDs)
    pool.map(get_illumina_dataframes, patientIDs)
