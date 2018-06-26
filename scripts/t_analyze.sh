import subprocess
import os
import shutil
targets=[[1106. 1142.]
 [1139. 1132.]
 [1275. 1180.]
 [1083. 1098.]
 [1275. 1158.]
 [1084. 1160.]
 [1204.  918.]
 [ 987. 1051.]
 [1128.  991.]
 [ 965. 1351.]]
for i in range(len(targets)):
try:
  os.mkdir('{}_{}'.format(targets[i][0],targets[i][1]))
except:
  pass
current='{}_{}'.format(targets[i][0],targets[i][1])
with open('run.py', 'rt') as fin:
  with open('to_qsub.py', 'wt') as fout:
    for line in fin:
      fout.write(line.replace('target', '{}'.format(targets[i])))
shutil.copyfile('to_qsub.py', '{}_{}/to_qsub.py'.format(targets[i][0],targets[i][1]))
shutil.copyfile('python.pbs', '{}_{}/python.pbs'.format(targets[i][0],targets[i][1]))
os.chdir(current)
os.system('qsub python.pbs')
os.chdir('../')