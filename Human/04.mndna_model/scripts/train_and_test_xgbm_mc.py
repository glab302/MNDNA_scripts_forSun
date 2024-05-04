import joblib
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import LabelEncoder
from xgboost.sklearn import XGBClassifier
from trainIO import *


def top_k_accuracy_per_label(y_true, y_pred_proba, k=2):
    # Get the unique labels
    labels = np.unique(y_true)

    # Initialize an empty dictionary to store the top-k accuracy for each label
    label_acc_dict = {}

    for label in labels:
        # Filter out the specific label
        indices = (y_true == label)
        y_true_label = y_true[indices]
        y_pred_proba_label = y_pred_proba[indices]

        # Get the top k predictions
        top_k_preds = np.argsort(y_pred_proba_label, axis=1)[:, -k:]

        # Check if the true labels are in top k predictions
        match_array = np.any(top_k_preds == y_true_label[:, None], axis=1)

        # Calculate top k accuracy
        top_k_acc = np.mean(match_array)

        # Add the top-k accuracy to the dictionary
        label_acc_dict[label] = top_k_acc

    return label_acc_dict


def calculate_score(y_true, y_pred):
    # 计算每个标签的准确率
    label_acc_dict = top_k_accuracy_per_label(y_true, y_pred, k=1)

    # 找到标签为1的准确率
    acc_label_1 = label_acc_dict[1] if 1 in label_acc_dict else 0

    # 计算所有标签的平均准确率
    avg_acc = sum(label_acc_dict.values()) / len(label_acc_dict)

    # 计算得分
    score = avg_acc

    return score, avg_acc, acc_label_1, label_acc_dict


def print_cm(cm, le):
    print("Confusion Matrix:")
    class_labels = [i for i in range(cm.shape[0])]
    print(',' + ','.join(str(le.inverse_transform([label])[0]) for label in class_labels))
    for encoded_label, row in enumerate(cm):
        original_label = le.inverse_transform([encoded_label])[0]
        # 计算真正例（True Positives）
        tp = cm[encoded_label, encoded_label]
        # 计算预测正例（Predicted Positives）
        pp = np.sum(cm[encoded_label, :])
        # 计算准确率（Accuracy）
        accuracy = tp / pp if pp > 0 else 0
        print(f"{original_label}," + ','.join(str(value) for value in row) + f",{accuracy:.2f}")


if __name__ == "__main__":
    print('XGBM Start: ' + str(datetime.datetime.now()))

    train_x, train_y, train_name, sc = loadfile('./learn.vector.train')
    test_x, test_y, test_name, sc = loadfile('./learn.vector.test', sc)
    le = LabelEncoder()

    all_y = np.concatenate([train_y, test_y])
    le.fit(all_y)

    train_y = le.transform(train_y)
    test_y = le.transform(test_y)

    # xgboost model
    xgbm = XGBClassifier(
        objective='multi:softprob',
        num_class=len(set(train_y)),
        learning_rate=0.1, n_estimators=350, max_depth=3, min_child_weight=1,
        seed=1234, subsample=0.8, colsample_bytree=1, gamma=0.7, reg_alpha=1, reg_lambda=1, n_jobs=12)

    fit2 = xgbm.fit(train_x, train_y)
    xgbm_predict = fit2.predict_proba(test_x)
    fit2.save_model('./xgbm.model.json')
    joblib.dump(sc, './xgbm.sc.pkl')

    xgbm_predict_train = fit2.predict(train_x)
    xgbm_predict_train_label = le.inverse_transform(xgbm_predict_train)
    xgbm_predict_prob_train = fit2.predict_proba(train_x)

    y_pred = fit2.predict(test_x)

    cm = confusion_matrix(test_y, y_pred)
    train_cm = confusion_matrix(train_y, xgbm_predict_train)

    print(top_k_accuracy_per_label(test_y, xgbm_predict, k=2))

    print_cm(cm, le)
    print_cm(train_cm, le)

    y_pred_label = le.inverse_transform(y_pred)
    y_test_label = le.inverse_transform(test_y)
    with open('./xgbm.predict_result', 'w') as file:
        for i in range(len(test_y)):
            file.write(f"{test_name[i]},{y_pred_label[i]},{y_test_label[i]}\n")
    # output each test sample's pred prob to file './xgbm.predict_prob_result'
    with open('./xgbm.predict_prob_result', 'w') as file:
        for i in range(len(xgbm_predict)):
            file.write(f"{test_name[i]},{','.join(str(value) for value in xgbm_predict[i])},{y_pred_label[i]}\n")
    # output each train sample's pred prob to file './xgbm.predict_prob_train'
    with open('./xgbm.predict_prob_train', 'w') as file:
        for i in range(len(xgbm_predict_prob_train)):
            file.write(f"{train_name[i]},"
                       + f"{','.join(str(value) for value in xgbm_predict_prob_train[i])},{xgbm_predict_train_label[i]}\n")
