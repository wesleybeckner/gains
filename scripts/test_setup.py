import subprocess
import os
import shutil

targets = [[813, 1310],
[434, 1360],
[334, 1342],
[397, 1274],
[741, 1021],
[713, 1024],
[525, 1444],
[434, 1396]]

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
      fout.write(line.replace('simplex_id', '{}'.format(i)))
shutil.copyfile('to_qsub.py', '{}_{}/to_qsub.py'.format(targets[i][0],targets[i][1]))
shutil.copyfile('python.pbs', '{}_{}/python.pbs'.format(targets[i][0],targets[i][1]))
os.chdir(current)
os.system('qsub python.pbs')
os.chdir('../')