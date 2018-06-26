import subprocess
import os
import pyinotify
import fnmatch
import re

# The watch manager stores the watches and provides operations on watches
wm = pyinotify.WatchManager()
mask = pyinotify.IN_CREATE  # watched events


class EventHandler(pyinotify.ProcessEvent):
    def process_IN_CREATE(self, event):
        # grab path, parent directory, and filename
        file = event.pathname.split("/")[-1]
        path = event.pathname.rsplit("/", 1)[0]
        parent_dir = event.pathname.rsplit("/", 2)[1]

        # submit boxmaker
        pair_id = file[1:]
        is_pair = find_pairs(path, pair_id)
        if len(is_pair) == 2:
            name = "C" + pair_id.split(".")[0] + "_A" + pair_id.split(".")[0]
            subprocess.call(['./spFF.sh %s' % name], shell=True)

        # submit en_min
        elif file == "conf.gro" and parent_dir == "boxmaker":
            name = re.search('C\S\S_A\S\S', path).group(0)
            subprocess.call(['./spMin.sh %s' % name], shell=True)

        # submit npt eq
        elif file == "confout.gro" and parent_dir == "min":
            name = re.search('C\S\S_A\S\S', path).group(0)
            subprocess.call(['./spEq.sh %s' % name], shell=True)

        # submit npt pr
        elif file == "confout.gro" and parent_dir == "equilibrate":
            name = re.search('C\S\S_A\S\S', path).group(0)
            subprocess.call(['./spNPT.sh %s' % name], shell=True)

        # run analysis
        elif file == "confout.gro" and parent_dir == "npt":
            name = re.search('C\S\S_A\S\S', path).group(0)
            subprocess.call(['./spAnalysis.sh %s' % name], shell=True)


def find_pairs(base, pattern):
    """Return list of files matching pattern in base folder"""
    return [n for n in fnmatch.filter(os.listdir(base), "?" + pattern) if
            os.path.isfile(os.path.join(base, n))]


handler = EventHandler()
notifier = pyinotify.Notifier(wm, handler)
wdd = wm.add_watch('../', mask, rec=True)
notifier.loop()