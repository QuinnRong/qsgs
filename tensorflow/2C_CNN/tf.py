import tensorflow as tf
import numpy as np
import os, sys
import time
import math

def get_label(filename):
    '''
    input: txt file
    idx1 val1
    idx2 val2
    ...
    output: an array of idx(string), an array of val(float)
    '''
    index = []
    label = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            index.append(line.split()[0])
            label.append(float(line.split()[1]))
    return np.array(index), np.array(label)

def get_cross_section(dire, dist, struct):
    '''
    input:
    direction: x, y or z
    distance: coordinate in the direction
    struct: array of shape (dim_x, dim_y, dim_z)
    output:
    2 dimensional cross section
    '''
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
    '''
    input:
    path + index[i] is a folder containing one stucture
    dire[i] + dist[i] represent a cross section
    output:
    (len(index), dim_h, dim_w, len(dire))
    '''
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
        image_data = []
        for i in range(len(dire)):
            image_data.append(get_cross_section(dire[i], dist[i], struct))
        data.append(np.array(image_data).transpose(1, 2, 0))
    return np.array(data)

def get_trai(path, dire, dist):
    '''
    input:
    path constains structures and result.txt
    dire, dist are lists determining cross sections used as channels
    output:
    label: (n,)
    data: (n, h, w, c)
    '''
    time_start=time.time()
    index, label = get_label(path + "/result.txt")
    data = get_data(path, index, dire, dist)
    log_file = open("log.txt", "a")
    print("\nX_train:", data.shape, file = log_file)
    print("y_train:", label.shape, file = log_file)
    time_end=time.time()
    print("time:", time_end-time_start, file = log_file)
    log_file.close()
    return data, label

def get_test(path, dire, dist):
    '''
    input:
    path constains structures and result.txt
    dire, dist are lists determining cross sections used as channels
    output:
    label: (n,)
    data: (n, h, w, c)
    '''
    time_start=time.time()
    index, label = get_label(path + "/result.txt")
    data = get_data(path, index, dire, dist)
    log_file = open("log.txt", "a")
    print("\nX_test:", data.shape, file = log_file)
    print("y_test:", label.shape, file = log_file)
    time_end=time.time()
    print("time:", time_end-time_start, file = log_file)
    log_file.close()
    return data, label

def complex_model(X, is_training, channels):
    # define weight
    Wconv1 = tf.get_variable("Wconv1", shape=[3, 3, channels, 32])
    bconv1 = tf.get_variable("bconv1", shape=[32])
    Wconv2 = tf.get_variable("Wconv2", shape=[3, 3, 32, 64])
    bconv2 = tf.get_variable("bconv2", shape=[64])
    Wconv3 = tf.get_variable("Wconv3", shape=[3, 3, 64, 128])
    bconv3 = tf.get_variable("bconv3", shape=[128])
    # fc1
    Wfc1   = tf.get_variable("Wfc1", shape=[4*4*128,256])
    # fc2
    Wfc2   = tf.get_variable("Wfc2", shape=[256,1])
    
    #define graph
    conv1 = tf.nn.conv2d(X,Wconv1,[1,2,2,1],padding="SAME") + bconv1
    conv1_norm = tf.layers.batch_normalization(conv1, training=is_training)
    ## 50->25
    relu1 = tf.nn.relu(conv1_norm)

    pool1 = tf.nn.max_pool(relu1,[1,2,2,1],padding="SAME",strides=[1,2,2,1])
    ## 25->13
    
    conv2 = tf.nn.conv2d(pool1,Wconv2,[1,2,2,1],padding="SAME") + bconv2
    conv2_norm = tf.layers.batch_normalization(conv2, training=is_training)
    ## 13->7
    relu2 = tf.nn.relu(conv2_norm)
    
    conv3 = tf.nn.conv2d(relu2,Wconv3,[1,2,2,1],padding="SAME") + bconv3
    conv3_norm = tf.layers.batch_normalization(conv3, training=is_training)
    ## 7->4
    relu3 = tf.nn.relu(conv3_norm)

    flat  = tf.reshape(relu3,[-1,4*4*128])
    fc1   = tf.matmul(flat,Wfc1)
    reluf = tf.nn.relu(fc1)
    fc2   = tf.matmul(reluf,Wfc2)
    
    return fc2[:, 0]

def run_test(sess, variables, feed_dict):
    loss, y_pred = sess.run(variables, feed_dict=feed_dict)
    np.savetxt("y_test_pred.txt", y_pred, fmt = '%.4f')
    return loss

def run_trai(sess, X, y, variables_trai, feed_dict_trai, variables_test, feed_dict_test, saver, epochs=1, batch_size=100):
    log_file = open("log.txt", "a")
    los_file = open("los.txt", "a")
    print("\nTraining", file = log_file)
    print("epoch batch     time    loss", file = log_file)
    # shuffle indicies
    Xd, yd = feed_dict_trai[X], feed_dict_trai[y]
    train_indicies = np.arange(Xd.shape[0])
    np.random.shuffle(train_indicies)
    np.savetxt("y_trai.txt", yd[train_indicies], fmt = '%.4f')
    # iterate over epoches
    time_start=time.time()
    for e in range(epochs):
        # keep track of losses and accuracy
        loss_trai = 0
        y_pred_trai = []
        for i in range(int(math.ceil(Xd.shape[0]/batch_size))):
            # generate indicies for the batch
            start_idx = (i*batch_size)%Xd.shape[0]
            idx = train_indicies[start_idx:start_idx+batch_size]
            # create a feed dictionary for this batch
            feed_dict_trai[X], feed_dict_trai[y] = Xd[idx,:], yd[idx]
            # get batch size
            actual_batch_size = yd[idx].shape[0]
            # have tensorflow compute losses and predictions and do gradient decent
            loss, y_pred, _ = sess.run(variables_trai,feed_dict=feed_dict_trai)
            loss_trai += loss*loss*actual_batch_size
            y_pred_trai.extend(list(y_pred))
            time_end=time.time()
            print("{0:5d} {1:5d} {2:8.1f} {3:7.3f}".format(e, i, time_end-time_start, loss), file = log_file)
        loss_trai = np.sqrt(loss_trai/Xd.shape[0])
        # testing results
        loss_test, y_pred_test = sess.run(variables_test, feed_dict=feed_dict_test)
        # save training and testing results
        print("Epoch {0:5d} training loss = {1:7.3f} testing loss = {2:7.3f}".format(e + 1, loss_trai, loss_test), file = los_file)
        np.savetxt("y_trai_pred.txt", y_pred_trai, fmt = '%.4f')
        np.savetxt("y_test_pred.txt", y_pred_test, fmt = '%.4f')
        # save the model
        if (e % 10 == 0):
            saver.save(sess, "Model/model.ckpt"+str(e))
            log_file.close()
            los_file.close()
            log_file = open("log.txt", "a")
            los_file = open("los.txt", "a")
    log_file.close()
    los_file.close()

def get_all_data():
    '''
    get training data and testing data
    '''
    dire = ["x", "y"]
    dist = [0, 0]
    test_data, test_label = get_test("../../structure/rand-50-0.250000_batch0", dire, dist)
    trai_data, trai_label = get_trai("../../structure/rand-50-0.250000_batch1", dire, dist)
    # normalize to mean 0, std 1
    mean, std = np.mean(trai_label), np.std(trai_label)
    trai_label = (trai_label - mean) / std
    test_label = (test_label - mean) / std
    # save results
    log_file = open("log.txt", "a")
    print("\nmean = {0:.3g}, std = {1:.3g}".format(mean, std), file = log_file)
    log_file.close()
    np.savetxt("y_test.txt", test_label, fmt = '%.4f')
    np.savetxt("y_trai.txt", trai_label, fmt = '%.4f')
    return test_data, test_label, trai_data, trai_label

def main():
    # testing data and training data
    test_data, test_label, trai_data, trai_label = get_all_data()
    # placeholders
    X = tf.placeholder(tf.float32, [None, 50, 50, test_data.shape[-1]])
    y = tf.placeholder(tf.float32, [None])
    is_training = tf.placeholder(tf.bool)
    # graph
    y_pred = complex_model(X, is_training, test_data.shape[-1])
    # loss
    loss = tf.sqrt(tf.losses.mean_squared_error(y_pred, y))
    # optimizer
    optimizer = tf.train.RMSPropOptimizer(1e-3)
    update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(update_ops):
        training = optimizer.minimize(loss)
    # saver
    saver = tf.train.Saver()
    # run tensorflow
    variables_test = [loss, y_pred]
    feed_dict_test = {X: test_data, y: test_label, is_training: False}
    variables_trai = [loss, y_pred, training]
    feed_dict_trai = {X: trai_data, y: trai_label, is_training: True}
    with tf.Session() as sess:
        sess.run(tf.global_variables_initializer())
        # restore results if needed
        # saver.restore(sess, "../run_8/Model/model.ckpt" + str(50))
        # check the initial loss
        log_file = open("log.txt", "a")
        print("\ninitial loss:", run_test(sess, variables_test, feed_dict_test), file = log_file)
        log_file.close()
        # train the model
        run_trai(sess, X, y, variables_trai, feed_dict_trai, variables_test, feed_dict_test, saver, 500, 100)

if __name__ == '__main__':
    sys.exit(main())
