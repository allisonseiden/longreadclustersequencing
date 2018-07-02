from PhasedData import PhasedData

patientIDs = ["1-00801", "1-01019", "1-03897", "1-04190", "1-04389", "1-04460",
                "1-04537", "1-05443", "1-05673", "1-05846"];
phased_data_objects = [];

for ID in patientIDs:
    phased_data_objects.append(PhasedData(ID));

for data_object in phased_data_objects:
    print(data_object);
