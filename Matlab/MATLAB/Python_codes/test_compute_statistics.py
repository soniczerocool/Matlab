'''test_compute_statistics'''


import visualise
import numpy as np
import SimpleITK as itk

file_path_seg = '/media/Expansion Drive/User/Research/Results/BRATS-2_KNN_ENDAUGUST/HG/0001/segmented_HG_0001.mha'
file_path_truth = '/media/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/HG/0001/VSD.Brain_3more.XX.XX.OT/VSD.Brain_3more.XX.XX.OT.6560.mha'
file_path_T1C = '/media/Expansion Drive/User/Research/Brain_data/BRATS-2/Image_Data/HG/0001/VSD.Brain.XX.O.MR_T1c/VSD.Brain.XX.O.MR_T1c.686.mha'


T1C = visualise.load_mha(file_path_T1C)