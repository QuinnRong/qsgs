import tensorflow as tf
# from skimage import io
import numpy as np
import os
import time
import math

def get_label(filename):
	index = []
	label = []
	with open(filename) as f:
		lines = f.readlines()
		for line in lines:
			index.append(line.split()[0])
			label.append(float(line.split()[1]))
	return np.array(index), np.array(label)

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

def get_data(path, index, dire, dist):
	data = []
	for id in index:
		newpath = os.path.join(path, id)
		struct = []
		count = 0
		for _ in os.listdir(newpath):
			count += 1
		for i in range(count):
			dat = "3D_" + id + "_" + str(i) + ".dat"
			temp = np.loadtxt(os.path.join(newpath, dat), dtype=int)
			struct.append(temp)
		struct = np.array(struct)
		image_data = get_cross_section(dire, dist, struct)
		data.append(image_data)
	return np.array(data)[:, :, :, np.newaxis]

def get_training(path, dire, dist):
	time_start=time.time()
	index, label = get_label(path + "/result.txt")
	data = get_data(path, index, dire, dist)
	log_file = open("log.txt", "a")
	print("X_train: ", data.shape, file = log_file)
	print("y_train: ", label.shape, file = log_file)
	time_end=time.time()
	print(time_end-time_start, file = log_file)
	log_file.close()
	return data, label

def get_testing(path, dire, dist):
	time_start=time.time()
	index, label = get_label(path + "/result.txt")
	data = get_data(path, index, dire, dist)
	log_file = open("log.txt", "a")
	print("X_test: ", data.shape, file = log_file)
	print("y_test: ", label.shape, file = log_file)
	time_end=time.time()
	print(time_end-time_start, file = log_file)
	log_file.close()
	return data, label

def complex_model(X):
	# define weight
	Wconv1 = tf.get_variable("Wconv1", shape=[3, 3, 1, 16])
	bconv1 = tf.get_variable("bconv1", shape=[16])
	Wconv2 = tf.get_variable("Wconv2", shape=[3, 3, 16, 32])
	bconv2 = tf.get_variable("bconv2", shape=[32])
	Wconv3 = tf.get_variable("Wconv3", shape=[3, 3, 32, 64])
	bconv3 = tf.get_variable("bconv3", shape=[64])

	# fc2
	Wfc2   = tf.get_variable("Wfc2",shape=[4*4*64,64])
	# fc3
	Wfc3   = tf.get_variable("Wfc3",[64,1])
	
	#define graph
	conv1 = tf.nn.conv2d(X,Wconv1,[1,2,2,1],padding="SAME") + bconv1
	## 50->25
	relu1 = tf.nn.relu(conv1)
	pool1 = tf.nn.max_pool(relu1,[1,2,2,1],padding="SAME",strides=[1,2,2,1])
	## 25->13
	
	conv2 = tf.nn.conv2d(pool1,Wconv2,[1,2,2,1],padding="SAME") + bconv2
	## 13->7
	relu2 = tf.nn.relu(conv2)
	
	conv3 = tf.nn.conv2d(relu2,Wconv3,[1,2,2,1],padding="SAME") + bconv3
	## 7->4
	relu3 = tf.nn.relu(conv3)

	flat   = tf.reshape(relu3,[-1,4*4*64])
	fc2   = tf.matmul(flat,Wfc2)
	reluf = tf.nn.relu(fc2)
	fc3   = tf.matmul(reluf,Wfc3)
	
	return fc3[:, 0]

def check(loss, y_out, data, label):
	with tf.Session() as sess:
		sess.run(tf.global_variables_initializer())
		smse, y_pred = sess.run([loss, y_out], feed_dict={X: data, y: label})
		log_file = open("log.txt", "a")
		print(np.sqrt(np.average(np.square(label - y_pred))), file = log_file)
		print(smse, file = log_file)
		log_file.close()

def run_model(session, loss_val, predict, Xd, yd, Xt, yt, epochs=1, batch_size=100, training=None):
	# shuffle indicies
	train_indicies = np.arange(Xd.shape[0])
	np.random.shuffle(train_indicies)
	np.savetxt("y_train.txt", yd[train_indicies], fmt = '%.4f')

	# setting up variables we want to compute (and optimizing)
	# if we have a training function, add that to things we compute
	variables = [loss_val, predict, predict]
	training_now = training is not None
	if training_now:
		variables[-1] = training

	time_start=time.time()
	for e in range(epochs):
		# keep track of losses and accuracy
		losses = []
		correct = []
		for i in range(int(math.ceil(Xd.shape[0]/batch_size))):
			# generate indicies for the batch
			start_idx = (i*batch_size)%Xd.shape[0]
			idx = train_indicies[start_idx:start_idx+batch_size]
			# create a feed dictionary for this batch
			feed_dict = {X: Xd[idx,:], y: yd[idx], is_training: training_now}
			# get batch size
			actual_batch_size = yd[idx].shape[0]

			# have tensorflow compute loss and correct predictions
			# and (if given) perform a training step
			loss, corr, _ = session.run(variables,feed_dict=feed_dict)
			losses.append(loss*loss*actual_batch_size)
			correct.extend(list(corr))
			print(len(correct))

			log_file = open("log.txt", "a")
			time_end=time.time()
			print(e, i, loss, time_end-time_start, file = log_file)
			log_file.close()

		np.savetxt("y_train_pred.txt", correct, fmt = '%.4f')
		loss_test, corr_test = session.run([loss_val, predict], feed_dict={X: Xt, y: yt})
		np.savetxt("y_test_pred.txt", corr_test, fmt = '%.4f')

		log_file = open("log.txt", "a")
		total_loss = np.sqrt(np.sum(losses)/Xd.shape[0])
		print("Epoch {0}, training loss = {1:.3g}, testing loss = {2:.3g}".format(e+1,total_loss,loss_test), file = log_file)
		log_file.close()

		if (e % 10 == 0):
			saver.save(sess, "Model/model.ckpt"+str(e))


test_data, test_label = get_testing("../../structure/rand-50-0.250000_batch0", "z", 0)
train_data, train_label = get_training("../../structure/rand-50-0.250000_batch1", "z", 0)

mean = np.mean(train_label)
std = np.std(train_label)
log_file = open("log.txt", "a")
print("mean = {0:.3g}, std = {1:.3g}".format(mean, std), file = log_file)
log_file.close()
train_label = (train_label - mean) / std
test_label = (test_label - mean) / std

np.savetxt("y_test.txt", test_label, fmt = '%.4f')
np.savetxt("y_train.txt", train_label, fmt = '%.4f')

X = tf.placeholder(tf.float32, [None, 50, 50, 1])
y = tf.placeholder(tf.float32, [None])
is_training = tf.placeholder(tf.bool)
y_out = complex_model(X)
loss = tf.sqrt(tf.losses.mean_squared_error(y_out, y))

check(loss, y_out, test_data, test_label)

optimizer = tf.train.AdamOptimizer(5e-5)
train_step = optimizer.minimize(loss)

saver = tf.train.Saver()

with tf.Session() as sess:
	sess.run(tf.global_variables_initializer())
	# saver.restore(sess, "../run_8/Model/model.ckpt" + str(50))

	log_file = open("log.txt", "a")
	print("Training", file = log_file)
	log_file.close()
	run_model(sess, loss, y_out, train_data, train_label, test_data, test_label, 200, 100, train_step)
