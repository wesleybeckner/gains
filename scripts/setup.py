import subprocess
import os
import shutil
targets = [[ 377, 1162],
       [ 811, 1238],
       [ 947, 1157],
       [ 436, 1202],
       [ 873,  995],
       [ 627, 1390],
       [ 889, 1408],
       [ 769, 1232],
       [ 975, 1024],
       [ 499, 1272]]

for i in range(len(targets)):
  try:
    os.mkdir("{}_{}".format(targets[i][0],targets[i][1]))
  except:
    pass
  current="{}_{}".format(targets[i][0],targets[i][1])
  with open("test_model.py", "rt") as fin:
    with open("to_qsub.py", "wt") as fout:
      for line in fin:
        fout.write(line.replace('target', '{}'.format(targets[i])))
  shutil.copyfile("to_qsub.py", "{}_{}/to_qsub.py".format(targets[i][0],targets[i][1]))
  shutil.copyfile("python.pbs", "{}_{}/python.pbs".format(targets[i][0],targets[i][1]))
  os.chdir(current)
  os.system("qsub python.pbs")
  os.chdir("../")

