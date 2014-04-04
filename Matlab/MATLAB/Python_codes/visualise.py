
import numpy as np 
import re
#import mha 
import os
import SimpleITK as itk

def load_mha(fn):
	I = itk.ReadImage(fn)
	height , width, depth = I.GetSize()
	Image = np.zeros((height, width, depth))

	for x in range(width):
		for y in range(height):
			for z in range(depth):
				Image[y,x,z] = I.GetPixel(y,x,z)

	
	return Image			



'''
def load(fn , size):
	""" loads the modalities of the given brain from textfiles"""

	height = size[0]
	width = size[1]
	depth = size[2]

	MASK = np.zeros((height, width, depth), dtype = np.int16)

	stream = open(os.path.expanduser(fn))
	lines = stream.readlines()	
	for l in lines:
		line = l.split()
		sColon = re.compile('[:]')
		ROW = int(round(float(sColon.split(line[4])[1]) * height))
		COL = int(round(float(sColon.split(line[5])[1]) * width))
		DEP = int(round(float(sColon.split(line[6])[1]) * depth))
		
		label = int(line[0])
		if label == 0:
			label = 10

		MASK[ROW,COL,DEP] = int(label)
	return MASK



img=mha.new()
img.read_mha('VSD.Brain.XX.O.MR_T1c.686_NOTCOMPRESSED.mha')	

MASK = load('interaction_HG_0001.txt', img.size)		

msk=mha.new(data=MASK, size=MASK.shape, spacing=img.spacing, offset=img.offset, data_type=img.data_type, direction_cosines=img.direction_cosines)
msk.write_mha('file_name.mha')
'''