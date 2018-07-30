from skimage import io
import shutil
import numpy as np
import os

def load_data(filename):
	data = np.loadtxt(filename)
	return data


def save_image(filename, data):
	io.imsave(filename, data)


def get_cross_section(dire, dist, struct):
	if dire == "x":
		return struct[:,dist,:]
	elif dire == "y":
		return struct[:,:,dist]
	elif dire == "z":
		return struct[dist,:,:]
	else:
		print("direction not correct!")
		exit()


def get_data(path, idx, dire, dist):
	newpath = os.path.join(path, str(idx))
	struct = []
	count = 0
	for _ in os.listdir(newpath):
		count += 1
	for i in range(count):
		dat = "3D_" + str(idx) + "_" + str(i) + ".dat"
		temp = load_data(os.path.join(newpath, dat))
		struct.append(temp)
	struct = np.array(struct)
	image_name = os.path.join(path, "image", str(idx) + "_" + dire + "_" + str(dist) + ".jpg")
	image_data = get_cross_section(dire, dist, struct)
	save_image(image_name, image_data)


def main(src, num, dim):
	if os.path.exists(os.path.join(src, "image")):
		shutil.rmtree(os.path.join(src, "image"))
	os.mkdir(os.path.join(src, "image"))
	for i in range(num):
		for j in range(dim):
			for d in ["x", "y", "z"]:
				get_data(src, i, d, j)


main("../structure/rand-50-0.250000_batch0", 10, 10)