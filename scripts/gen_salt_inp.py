#script to generate salt.inp file
import re
log = pd.read_csv("salt_log.csv")
salt = []
for i in range(log.shape[0]):
    prediction = re.findall("\d+\.\d+", log["Model Prediction"][i])
    if i < 10:
        salt.append("C0{},A0{},{:.3}".format(i + 1, i + 1,
                                             float(prediction[0]) /
                                             1e3 - 0.2))
    else:
        salt.append("C{},A{},{:.3}".format(i + 1, i + 1,
                                           float(prediction[0]) /
                                             1e3 - 0.2))
with open('salt.inp', mode='wt', encoding='utf-8') as myfile:
    myfile.write('\n'.join(salt))
