from skimage import io
import shutil
import numpy as np
import os

def load_data(filename):
	data = np.loadtxt(filename)
	return data


def save_image(filename, data):
	io.imsave(filename, data)


def main(src, num, dim):
	if os.path.exists(os.path.join(src, "image")):
		shutil.rmtree(os.path.join(src, "image"))
	os.mkdir(os.path.join(src, "image"))
	for i in range(num):
		for j in range(dim):
			name_data = os.path.join(src, str(i), "3D_" + str(i) + "_" + str(j) + ".dat")
			data = load_data(name_data)
			name_image = os.path.join(src, "image", str(i) + "_" + str(j) + ".jpg")
			save_image(name_image, data)

# main("structure/core-50-0.250000", 1, 10)
main("structure/parallel-20-0.200000", 1, 20)
main("structure/serial-20-0.200000", 1, 20)