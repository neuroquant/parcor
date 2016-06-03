from pymatbridge import Matlab
mlab = Matlab()
mlab.start()
res = mlab.run_code('clear; main.m');
mlab.stop()
