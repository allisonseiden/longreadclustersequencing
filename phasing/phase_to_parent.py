from PhasedData import PhasedData;
import pandas as pd;
import numpy as np;
import multiprocessing as mp;



# patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389"];
patientIDs = ['1-04460', '1-04537', '1-05443', '1-05673', '1-05846'];

def run_phase_to_parent(ID):
    patient = PhasedData(ID);
    patient.create_vcf_dictionary();
    patient.create_dnvs_dictionary();
    patient.fill_bounds_dictionary();
    patient.find_variants_for_phasing(7);
    patient.assign_to_parent();
    patient.convert_to_dataframe();

if __name__ == '__main__':
    pool = mp.Pool(processes=5);
    pool.map(run_phase_to_parent, patientIDs);
