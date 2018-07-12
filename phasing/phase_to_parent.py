from PhasedData import PhasedData;
import pandas as pd;
import numpy as np;
import multiprocessing as mp;



patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389"];

def run_phase_to_parent(ID):
    patient = PhasedData(ID);
    patient.create_vcf_dictionary();
    patient.create_dnvs_dictionary();
    patient.fill_bounds_dictionary();
    patient.find_variants_for_phasing();
    patient.assign_to_parent();
    patient.convert_to_dataframe();

if __name__ == '__main__':
    pool = mp.Pool(processes=5);
    pool.map(run_phase_to_parent, patientIDs);
