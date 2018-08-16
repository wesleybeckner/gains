import subprocess
import os
import shutil

targets = [[1570, 899],
[306, 1407],
[233, 1199],
[413, 1099],
[1325, 931],
[504, 1034],
[584, 1532],
[346, 1471]]

for i in range(len(targets)):
  try:
    os.mkdir('{}_{}'.format(targets[i][0],targets[i][1]))
  except:
    pass
  current='{}_{}'.format(targets[i][0],targets[i][1])
  with open('run.py', 'rt') as fin:
    with open('run2.py', 'wt') as fout:
      for line in fin:
        fout.write(line.replace('target', '{}'.format(targets[i])))
  with open('run2.py', 'rt') as fin:
    with open('to_qsub.py', 'wt') as fout:
      for line in fin:
        fout.write(line.replace('simplex_id', '{}'.format(i)))
  shutil.copyfile('to_qsub.py', '{}_{}/to_qsub.py'.format(targets[i][0],targets[i][1]))
  shutil.copyfile('python.pbs', '{}_{}/python.pbs'.format(targets[i][0],targets[i][1]))
  os.chdir(current)
  os.system('qsub python.pbs')
  os.chdir('../')