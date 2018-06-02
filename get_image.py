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

# main("50-parallel", 1, 50)
# main("20-0.200000-serial", 1, 20)
# main("20-0.400000-0.020000", 1, 20)
main("100-0.250000-0.010000", 1, 100)
# main("100-0.400000-0.020000", 1, 100)
