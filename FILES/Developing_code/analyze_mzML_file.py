import json

import pandas as pd
import os


from pyopenms import *#library used for processing of MS 

from massql import msql_fileloading #used to load metabolite data

from comparison.find_triplet import find_triplet #delete I dont need to find triplets

def _determine_scan_polarity_pyteomics_mzML(spec):
    """
    Gets an enum for positive and negative polarity, for pyteomics

    Args:
        spec ([type]): [description]

    Returns:
        [type]: [description]
    """
    polarity = 0

    if "negative scan" in spec:
        polarity = 2
    if "positive scan" in spec:
        polarity = 1

    return polarity

def feature_finding(input_filename):
    options = PeakFileOptions()
    options.setMSLevels([1])
    fh = MzMLFile()
    fh.setOptions(options)

    # Load data
    input_map = MSExperiment()
    fh.load(input_filename, input_map)
    input_map.updateRanges()

    ff = FeatureFinder()
    ff.setLogType(LogType.CMD)

    # Run the feature finder
    name = "centroided"
    features = FeatureMap()
    seeds = FeatureMap()
    params = FeatureFinder().getParameters(name)
    ff.run(name, input_map, features, params, seeds)

    features.setUniqueIds()
    fh = FeatureXMLFile()

    output_filename = input_filename.replace(".mzMl", "_output.featureXML")
    fh.store(output_filename, features)
    return features


def triplet_extraction(input_filename):
    ms1_df, ms2_df = msql_fileloading.load_data(input_filename, cache='feather')

    featurexml_filename = input_filename.replace(".mzMl", "_output.featureXML")

    current_directory = os.getcwd()

    featurexml_file = os.path.join(current_directory, featurexml_filename)

    if os.path.exists(featurexml_file):
        features = FeatureMap()
        FeatureXMLFile().load(featurexml_file, features)
    else:
        features = feature_finding(input_filename)



    spectra = []
    i = 0
    duplas = []
    maxMs2i = []
    features_scans = []
    for f in features:
        print(ms2_df.loc[(
                                            (f.getMZ() - 0.1 < ms2_df['precmz']) & (ms2_df['precmz'] < f.getMZ() + 0.1)) & (
                                            (f.getRT() / 60 - 0.1 < ms2_df['rt']) & (ms2_df['rt'] < f.getRT() / 60 + 0.1))][
                             'scan'].unique())
        matched_scans = (ms2_df.loc[(
                                            (f.getMZ() - 0.1 < ms2_df['precmz']) & (ms2_df['precmz'] < f.getMZ() + 0.1)) & (
                                            (f.getRT() / 60 - 0.1 < ms2_df['rt']) & (ms2_df['rt'] < f.getRT() / 60 + 0.1))][
                             'scan'].unique())
        features_scans.append(matched_scans)
        for m in matched_scans:
            maxMs2i.append(ms2_df.loc[ms2_df['scan'] == m]['i'].max())
        Scan_intensity_dict = {
            'scan': matched_scans,
            'maxMs2i': maxMs2i
        }
        Scan_intensity_df = pd.DataFrame(Scan_intensity_dict)
        Scan_intensity_df_sorted = Scan_intensity_df.sort_values(by='maxMs2i')
        scans = Scan_intensity_df_sorted['scan']
        if len(scans) > 1:
            duplas.append(scans[-2:])

        spectra.clear()
        maxMs2i.clear()
        i = i + 1

    triplets = []
    i = len(duplas)
    for dupla in duplas:
        print(i)
        triplet, comparison_scores = find_triplet(dupla, features_scans, ms2_df)
        print(triplet)
        if(triplet!=[]):
            triplets_dict = {
                'dupla': dupla,
                'triplet': triplet,
                'scores' : comparison_scores
            }
            triplets.append(triplets_dict)
        i = i - 1

    return triplets

def main():
    # Ruta del archivo mzML que deseas analizar
    mzml_file_path = r"C:\Users\maipa\.spyder-py3\160528_FHS_Eicosanoid_Plate_13_03.mzML"
    # Ejecutar la extracción de tripletes
    print("Ejecutando análisis de metabolitos en:", mzml_file_path)
    triplets = triplet_extraction(mzml_file_path)

    # Mostrar los resultados obtenidos
    if triplets:
        print(f"Se encontraron {len(triplets)} tripletes.")
        for t in triplets:
            print(t)
    else:
        print("No se encontraron tripletes.")

if __name__ == "__main__":
    main()
    