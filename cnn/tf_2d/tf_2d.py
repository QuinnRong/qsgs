import tensorflow as tf
import numpy as np
import os, sys
import time
import math

# global parameters
resolution = 100    # input resulation
dp_parm = 1.0       # keep probability of dropout

class Param:
    def __init__(self, label, struc, start, end):
        self.label = label
        self.struc = struc
        self.start = start
        self.end   = end

valid_param = Param("../../../fenics/run_1_valid/output", "../../../utility/run_1_valid/output", 0, 199)
train_param = Param("../../../fenics/run_2_train/output", "../../../utility/run_2_train/output", 0, 799)

def get_label(filename, start, end):
    '''
    input: txt file, range from start to end
    idx1 val1
    idx2 val2
    ...
    output: an array of idx(string), an array of val(float)
    '''
    index = []
    label = []
    with open(filename) as f:
        lines = f.readlines()
        if len(lines) <= end:
            sys.exit("index out of range!")
        for i in range(start, end + 1):
            index.append(lines[i].split()[0])
            label.append(float(lines[i].split()[1]))
    return np.array(index), np.array(label)

def get_cross_section(dire, dist, struct):
    '''
    input:
    direction: x, y or z
    distance: coordinate in the direction
    struct: array of shape (dim_z, dim_y, dim_x)
    output:
    2 dimensional cross section
    '''
    if dire == "x":
        return struct[:,:,dist]
    elif dire == "y":
        return struct[:,dist,:]
    elif dire == "z":
        return struct[dist,:,:]
    else:
        sys.exit("direction not correct!")

def get_struc(path, index, dire, dist):
    '''
    input:
    path + index[i] defines a stucture file
    dire[i] + dist[i] represent a cross section
    output:
    (len(index), dim_h, dim_w, len(dire))
    '''
    struc = []
    for id in index:
        filename = "3D_" + id + ".dat"
        data = np.loadtxt(os.path.join(path, filename), dtype=int)
        dim_x = data.shape[1]
        if dim_x != resolution:
            sys.exit("resolution not match!")
        if data.shape[0] != dim_x**2:
            sys.exit("input is not a cubic!")
        data = data.reshape((resolution, resolution, resolution))
        image_data = []
        for i in range(len(dire)):
            image_data.append(get_cross_section(dire[i], dist[i], data))
        struc.append(np.array(image_data).transpose(1, 2, 0))
    return np.array(struc)

def get_tf_input(labelpath, structpath, start, end, dire, dist):
    '''
    input:
    path labelpath containing result.txt and structpath containing 3D_n.dat
    start and end defines index range
    output:
    label: (end-start+1,)
    struc: (end-start+1, dim, dim, len(dire))
    '''
    index, label = get_label(labelpath + "/result.txt", start, end)
    struc = get_struc(structpath, index, dire, dist)
    return struc, label

def get_all_data():
    '''
    get training data and validating data
    '''
    dire = ["x", "y", "x", "y"]
    dist = [0, 0, 50, 50]
    # validating data
    log_file = open("log.txt", "a")
    time_start=time.time()
    valid_struc, valid_label = get_tf_input(valid_param.label, valid_param.struc, valid_param.start, valid_param.end, dire, dist)
    time_end=time.time()
    print("\nX_valid: ", valid_struc.shape, file = log_file)
    print("y_valid: ", valid_label.shape, file = log_file)
    print(time_end-time_start, file = log_file)
    log_file.close()
    # training data
    log_file = open("log.txt", "a")
    time_start=time.time()
    train_struc, train_label = get_tf_input(train_param.label, train_param.struc, train_param.start, train_param.end, dire, dist)
    time_end=time.time()
    print("\nX_train: ", train_struc.shape, file = log_file)
    print("y_train: ", train_label.shape, file = log_file)
    print(time_end-time_start, file = log_file)
    log_file.close()
    # normalize to mean 0, std 1
    mean, std = np.mean(train_label), np.std(train_label)
    train_label = (train_label - mean) / std
    valid_label = (valid_label - mean) / std
    # save results
    log_file = open("log.txt", "a")
    print("mean = {0:.6g}, std = {1:.6g}".format(mean, std), file = log_file)
    log_file.close()
    return valid_struc, valid_label, train_struc, train_label, mean, std

def complex_model(X, is_training, channels, keep_prob, name):
    # define weight
    Wconv1 = tf.get_variable("Wconv1", shape=[5, 5, channels, 64])
    bconv1 = tf.get_variable("bconv1", shape=[64])
    Wconva = tf.get_variable("Wconva", shape=[3, 3, 64, 64])
    bconva = tf.get_variable("bconva", shape=[64])
    Wconv2 = tf.get_variable("Wconv2", shape=[3, 3, 64, 128])
    bconv2 = tf.get_variable("bconv2", shape=[128])
    Wconvb = tf.get_variable("Wconvb", shape=[3, 3, 128, 128])
    bconvb = tf.get_variable("bconvb", shape=[128])
    Wconv3 = tf.get_variable("Wconv3", shape=[3, 3, 128, 256])
    bconv3 = tf.get_variable("bconv3", shape=[256])
    # fc1
    Wfc1   = tf.get_variable("Wfc1", shape=[4*4*256,512])
    # fc2
    Wfc2   = tf.get_variable("Wfc2", shape=[512,1])
    
    #define graph
    conv1 = tf.nn.conv2d(X,Wconv1,[1,2,2,1],padding="SAME") + bconv1
    conv1_norm = tf.layers.batch_normalization(conv1, training=is_training)
    relu1 = tf.nn.relu(conv1_norm)
    ## 100->50
    conva = tf.nn.conv2d(relu1,Wconva,[1,1,1,1],padding="SAME") + bconva
    conva_norm = tf.layers.batch_normalization(conva, training=is_training)
    relua = tf.nn.relu(conva_norm)

    pool1 = tf.nn.max_pool(relua,[1,2,2,1],padding="SAME",strides=[1,2,2,1])
    ## 50->25

    conv2 = tf.nn.conv2d(pool1,Wconv2,[1,2,2,1],padding="SAME") + bconv2
    conv2_norm = tf.layers.batch_normalization(conv2, training=is_training)
    relu2 = tf.nn.relu(conv2_norm)
    ## 25->13
    convb = tf.nn.conv2d(relu2,Wconvb,[1,1,1,1],padding="SAME") + bconvb
    convb_norm = tf.layers.batch_normalization(convb, training=is_training)
    relub = tf.nn.relu(convb_norm)

    pool2 = tf.nn.max_pool(relub,[1,2,2,1],padding="SAME",strides=[1,2,2,1])
    ## 13->7
 
    conv3 = tf.nn.conv2d(pool2,Wconv3,[1,2,2,1],padding="SAME") + bconv3
    conv3_norm = tf.layers.batch_normalization(conv3, training=is_training)
    relu3 = tf.nn.relu(conv3_norm)
    ## 7->4

    dp   = tf.nn.dropout(relu3, keep_prob)
    flat  = tf.reshape(dp,[-1,4*4*256])
    fc1   = tf.matmul(flat,Wfc1)
    reluf = tf.nn.relu(fc1)
    fc2   = tf.matmul(reluf,Wfc2, name=name)
    
    return fc2[:, 0]

def run_valid(sess, variables, feed_dict):
    loss, y_pred = sess.run(variables, feed_dict=feed_dict)
    return loss

def run_train(sess, X, y, variables_train, feed_dict_train, variables_valid, feed_dict_valid, saver, mean, std, epochs, batch_size):
    log_file = open("log.txt", "a")
    los_file = open("los.txt", "a")
    print("\nTraining", file = log_file)
    print("epoch batch     time     loss", file = log_file)
    # shuffle indicies
    Xd, yd = feed_dict_train[X], feed_dict_train[y]
    train_indicies = np.arange(Xd.shape[0])
    np.random.shuffle(train_indicies)
    np.savetxt("y_train.txt", yd[train_indicies], fmt = '%.5f')
    np.savetxt("y_valid.txt", feed_dict_valid[y], fmt = '%.5f')
    # iterate over epoches
    time_start=time.time()
    # mae_min = 0.5
    for e in range(1, epochs+1):
        # keep track of losses and accuracy
        loss_train = 0
        y_pred_train = []
        for i in range(int(math.ceil(Xd.shape[0]/batch_size))):
            # generate indicies for the batch
            start_idx = (i*batch_size)%Xd.shape[0]
            idx = train_indicies[start_idx:start_idx+batch_size]
            # create a feed dictionary for this batch
            feed_dict_train[X], feed_dict_train[y] = Xd[idx,:], yd[idx]
            # get batch size
            actual_batch_size = yd[idx].shape[0]
            # have tensorflow compute losses and predictions and do gradient decent
            loss, y_pred, _ = sess.run(variables_train,feed_dict=feed_dict_train)
            loss_train += loss*loss*actual_batch_size
            y_pred_train.extend(list(y_pred))
            time_end=time.time()
            print("{0:5d} {1:5d} {2:8.1f} {3:7.3f}".format(e, i, time_end-time_start, loss), file = log_file)
            log_file.close()
            log_file = open("log.txt", "a")
        loss_train = np.sqrt(loss_train/Xd.shape[0])
        # validating results
        loss_valid, y_pred_valid = sess.run(variables_valid, feed_dict=feed_dict_valid)
        mae = np.mean(np.abs(y_pred_valid - feed_dict_valid[y]))
        mae = mae * std / mean
        # save training and validating results
        print("Epoch{0:5d} trainloss={1:6.3f} validloss={2:6.3f} mae={3:8.5f}".format(e,loss_train,loss_valid,mae), file=los_file)
        los_file.close()
        los_file = open("los.txt", "a")
        # save the model
        # if (e > 100 and mae < mae_min):
        if (e >= 400 and e % 20 == 0)
            # mae_min = mae
            saver.save(sess, "Model/model.ckpt"+str(e))
            np.savetxt("y_valid_"+str(e)+".txt", y_pred_valid, fmt = '%.5f')
            np.savetxt("y_train_"+str(e)+".txt", y_pred_train, fmt = '%.5f')
    log_file.close()
    los_file.close()

def main():
    # testing data and training data
    valid_struc, valid_label, train_struc, train_label, mean, std = get_all_data()
    # placeholders
    X = tf.placeholder(tf.float32, [None, resolution, resolution, valid_struc.shape[-1]], name="X")
    y = tf.placeholder(tf.float32, [None], name="y")
    is_training = tf.placeholder(tf.bool, name="is_training")
    keep_prob = tf.placeholder(tf.float32, name="keep_prob")
    # graph
    y_pred = complex_model(X, is_training, valid_struc.shape[-1], keep_prob, name="y_pred")
    # loss
    loss = tf.sqrt(tf.losses.mean_squared_error(y_pred, y), name="loss")
    # optimizer
    optimizer = tf.train.RMSPropOptimizer(1e-3)
    update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(update_ops):
        training = optimizer.minimize(loss)
    # saver
    saver = tf.train.Saver()
    # run tensorflow
    variables_valid = [loss, y_pred]
    feed_dict_valid = {X: valid_struc, y: valid_label, is_training: True, keep_prob: 1}
    variables_train = [loss, y_pred, training]
    feed_dict_train = {X: train_struc, y: train_label, is_training: True, keep_prob: dp_parm}
    with tf.Session() as sess:
        # restore weights or initialize weights
        # saver.restore(sess, "./Model/model.ckpt" + str(50))
        sess.run(tf.global_variables_initializer())
        # check the initial loss
        log_file = open("log.txt", "a")
        print("\nvalid loss:", run_valid(sess, variables_valid, feed_dict_valid), file = log_file)
        log_file.close()
        # train the model
        run_train(sess, X, y, variables_train, feed_dict_train, variables_valid, feed_dict_valid, saver, mean, std, 500, 100)

if __name__ == '__main__':
    sys.exit(main())
