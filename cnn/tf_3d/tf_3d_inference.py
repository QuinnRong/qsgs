import tensorflow as tf
import re, os, sys
import time
import numpy as np

# global parameters
resolution = 100
models = [100, 200, 300, 400, 500]
model_file = "./Model/model.ckpt"

class Param:
    def __init__(self, label, struc, start, end):
        self.label = label
        self.struc = struc
        self.start = start
        self.end   = end

test_param = Param("../../../fenics/run_3_test/output", "../../../utility/run_3_test/output", 0, 199)

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

def get_struc(path, index):
    '''
    input:
    path + index[i] defines a stucture file
    output:
    (len(index), dim_z, dim_y, dim_x, 1)
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
        struc.append(data)
    return np.array(struc)[:, :, :, :, np.newaxis]

def get_tf_input(labelpath, structpath, start, end):
    '''
    input:
    path labelpath containing result.txt and structpath containing 3D_n.dat
    start and end defines index range
    output:
    label: (end-start+1,)
    struc: (end-start+1, dim, dim, dim, 1)
    '''
    index, label = get_label(labelpath + "/result.txt", start, end)
    struc = get_struc(structpath, index)
    return struc, label

def get_mean_std(filename, nline):
    '''
    input:
    nth line in the file is, for example,
    mean = 2.5549, std = 0.871817
    output:
    float(2.5549), float(0.871817)
    '''
    with open(filename) as f:
        for i in range(nline - 1):
            f.readline()
        line = f.readline()
        linelist = re.split("[=,]", line)
    return float(linelist[1]), float(linelist[3])

def get_all_data():
    '''
    get testing data
    '''
    time_start=time.time()
    test_struc, test_label = get_tf_input(test_param.label, test_param.struc, test_param.start, test_param.end)
    time_end=time.time()
    log_file = open("log.txt", "a")
    print("\nX_test: ", test_struc.shape, file=log_file)
    print("y_test: ", test_label.shape, file=log_file)
    print(time_end-time_start, file=log_file)
    log_file.close()
    return test_struc, test_label

def main():
    # testing data
    test_struc, test_label = get_all_data()
    mean, std = get_mean_std("log.txt", 9)
    np.savetxt("y_test.txt", test_label, fmt = '%.5f')
    test_label = (test_label - mean) / std 
    # run tensorflow
    maes = []
    rmses = []
    for m in models:
        with tf.Session() as sess:
            saver = tf.train.import_meta_graph(model_file + str(m) + ".meta")
            saver.restore(sess, model_file + str(m))

        graph = tf.get_default_graph()
        loss = graph.get_operation_by_name("loss").outputs[0]                   # loss:0
        y_pred = graph.get_operation_by_name("y_pred").outputs[0]               # y_pred:0
        X = graph.get_operation_by_name("X").outputs[0]                         # X:0
        y = graph.get_operation_by_name("y").outputs[0]                         # y:0
        is_training = graph.get_operation_by_name("is_training").outputs[0]     # is_training:0

        variables_test = [loss, y_pred[:,0]]
        feed_dict_test = {X: test_struc, y: test_label, is_training: True}

            loss, y_pred = sess.run(variables_test, feed_dict=feed_dict_test)
            np.savetxt("y_test_"+str(m)+".txt", y_pred*std+mean, fmt = '%.5f')
            mae = np.mean(np.abs(y_pred - test_label)) * std / mean
            maes.append(mae)
            rmse = np.sqrt(np.mean(np.square((y_pred-test_label)*std/(test_label*std+mean))))
            rmses.append(rmse)

            log_file = open("log.txt", "a")
            print("\nloss = ", loss, file=log_file)
            print("mean absolute error = ", mae, file=log_file)
            print("root mean square error = ", rmse, file=log_file)
            log_file.close()

    log_file = open("log.txt", "a")
    print("\nmean absolute error: ", np.mean(maes), " ", np.std(maes), file=log_file)
    print("root mean square error: ", np.mean(rmses), " ", np.std(rmses), file=log_file)
    log_file.close()

if __name__ == '__main__':
    sys.exit(main())
