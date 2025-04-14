# -*- coding: utf-8 -*-
# main.py
# -*- coding: utf-8 -*-
"""
PARTE A
Aplico filtro a la señal y obtengo el baseline

PARTE B
1. Calculo para cada intensidad su frecuencia correspondiente f(i)
2. Creo una intensidad newI=Ii-f(i)xA, siiendo A la amplitud para cada mz
3. Me genero mi oscilación: le resto a cada intensidad de la señal original(filtrada) 
la frecuencia para esa intensidad x la amplitud en esa intensidad i

"""
from processing.processor import process_file
import os

if __name__ == "__main__":
    plot_count = 0  # Global variable for validations

    file_path = "C:/Users/maipa/repos/oscilations_TFG_repository/FILES/TESTS/MSCONVERT_CTRL_103_01_c_afterreboot_original.mzML"
    save_as = "C:/Users/maipa/repos/oscilations_TFG_repository/FILES/TESTS/CTRL_103_01_c_afterreboot_corrected.mzML"

    if os.path.exists(save_as):
        print(f"El archivo {save_as} ya existe. Eliminándolo para reemplazarlo...")
        os.remove(save_as)

    process_file(file_path, save_as, window_length=11, filter_order=5)


