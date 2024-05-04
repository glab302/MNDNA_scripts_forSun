import datetime, math, os, random
import numpy as np
from sklearn import preprocessing


def loadfile(path, scaler='none'):
    x = []
    y = []
    name = []
    with open(path, 'r') as f:
        n, m = [int(x) for x in f.readline().split(' ')]
        for i in range(n):
            line = f.readline()
            ss = line.split(' ')
            y.append(int(ss[0]))
            name.append(ss[1])
            x.append([float(t) for t in ss[2:]])
    #            a = [float(x) for x in line.split(' ')]
    #            y.append(a[0])
    #            x.append(a[2:])
    #            name.append(a[1])
    x1 = np.array(x)
    #    if scaler == 'none':
    return x1, np.array(y), name, scaler


#    scaler = preprocessing.StandardScaler().fit(x1)
#    return scaler.transform(x1), np.array(y), name, scaler

def get_val(X, Y):
    if Y == 0:
        return "N/A"
    return str(X * 100 / Y)


def write_array_to_file(arr, fil):
    if arr.ndim == 2:
        a2 = arr[0]
    elif arr.ndim == 1:
        a2 = arr
    else:
        return
    s = '['
    for C in a2:
        s += str(C) + ', '
    s += ']\n'
    fil.write(s)


def self_check(T, stand):
    l = len(stand)
    A, B, C, D = 0, 0, 0, 0
    test = []
    for i in range(l):
        if not (isinstance(T[i], np.ndarray)):
            T1 = T[i]
            T0 = 1 - T1
        elif len(T[i]) > 1:
            T0 = T[i][0]
            T1 = T[i][1]
        else:
            T1 = T[i][0]
            T0 = 1 - T1
        if T0 > T1:
            test.append(0)
        else:
            test.append(1)

    for i in range(l):
        if stand[i] == 0 and test[i] == 0:
            D += 1
        elif stand[i] == 1 and test[i] == 1:
            A += 1
        elif stand[i] == 0 and test[i] == 1:
            B += 1
        elif stand[i] == 1 and test[i] == 0:
            C += 1
    print('Training Data Sample Positive, Predict Positive: ' + str(A) + '\n')
    print('Training Data Sample Negative, Predict Positive: ' + str(B) + '\n')
    print('Training Data Sample Positive, Predict Negative: ' + str(C) + '\n')
    print('Training Data Sample Negative, Predict Negative: ' + str(D) + '\n')


def writefile(stand, name, T, result_path, stat_path):
    l = len(stand)
    A, B, C, D = 0, 0, 0, 0
    test = []
    pos = []
    neg = []

    with open(result_path, 'w') as result:
        result.write(str(l) + '\n')
        for i in range(l):
            if not (isinstance(T[i], np.ndarray)):
                T1 = T[i]
                T0 = 1 - T1
            elif len(T[i]) > 1:
                T0 = T[i][0]
                T1 = T[i][1]
            else:
                T1 = T[i][0]
                T0 = 1 - T1
            result.write(str(name[i]) + ' ' + str(T1) + ' ')
            if stand[i] == 1:
                pos.append(T1)
            else:
                neg.append(T1)
            if T0 > T1:
                test.append(0)
                result.write('0')
            else:
                test.append(1)
                result.write('1')
            result.write('\n')

    pos.sort()
    neg.sort()
    auc = 0
    for pos_p in pos:
        for neg_p in neg:
            if pos_p > neg_p:
                auc += 1.0
    auc = get_val(auc, len(pos) * len(neg))

    for i in range(l):
        if stand[i] == 0 and test[i] == 0:
            D += 1
        elif stand[i] == 1 and test[i] == 1:
            A += 1
        elif stand[i] == 0 and test[i] == 1:
            B += 1
        elif stand[i] == 1 and test[i] == 0:
            C += 1
    with open(stat_path, 'w') as f:
        f.write('AUC: ' + str(auc) + '\n')
        f.write('Sample Positive, Predict Positive: ' + str(A) + '\n')
        f.write('Sample Negative, Predict Positive: ' + str(B) + '\n')
        f.write('Sample Positive, Predict Negative: ' + str(C) + '\n')
        f.write('Sample Negative, Predict Negative: ' + str(D) + '\n')
        f.write('Accuracy: ' + get_val(A + D, A + B + C + D) + '%' + '\n')
        f.write('PPV: ' + get_val(A, A + B) + '%' + '\n')
        f.write('NPV: ' + get_val(D, C + D) + '%' + '\n')
        f.write('Sensitivity: ' + get_val(A, A + C) + '%' + '\n')
        f.write('Specificity: ' + get_val(D, B + D) + '%' + '\n')
        f.write('Precision: ' + get_val(A, A + B) + '%' + '\n')