from xgboost.sklearn import XGBClassifier
import joblib
from trainIO import *

if __name__ == "__main__":
    print('XGBM Start: ' + str(datetime.datetime.now()))

    train_x, train_y, train_name, sc = loadfile('./learn.vector.train')
    test_x, test_y, test_name, sc = loadfile('./learn.vector.test', sc)
    xgbm = XGBClassifier(learning_rate=0.1, n_estimators=190, max_depth=3, min_child_weight=1, seed=1234, subsample=1, colsample_bytree=1, gamma=0.7, reg_alpha=1, reg_lambda=1, n_jobs=12)
    fit2 = xgbm.fit(train_x, train_y)
    fit2.save_model('./xgbm.model.json')
    joblib.dump(sc, './xgbm.sc.pkl')

    xgbm_predict = fit2.predict_proba(test_x)
    xgbm_predict_train = fit2.predict_proba(train_x)
    self_check(fit2.predict_proba(train_x), train_y)
    writefile(test_y, test_name, xgbm_predict, './xgbm.predict_result', './xgbm.predict_stat')
    # writefile(train_y, train_name, xgbm_predict_train, './xgbm.predict_result_train', './xgbm.predict_stat_train')

    print('End: ' + str(datetime.datetime.now()))