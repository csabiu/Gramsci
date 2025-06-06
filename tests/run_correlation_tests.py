import os
import subprocess
import math
import numpy as np

def run(cmd):
    print(cmd)
    subprocess.check_call(cmd, shell=True)

def assert_close(val, expected, tol, msg):
    if abs(val - expected)/expected > tol:
        raise AssertionError(f"{msg}: got {val} expected {expected}")

def main():
    bindir = os.path.join('..', 'bin')
    gal = os.path.join('..', 'example', 'test.gal')
    ran = os.path.join('..', 'example', 'test.ran')

    out2 = 'tmp_2pcf.out'
    cmd2 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 10 -wgt -nmu 10 -out {out2} -2pcf"
    run(cmd2)

    tmp=np.loadtxt(out2,skiprows=1)

    DD=np.mean(tmp[:,4])
    RR=np.mean(tmp[:,5])

    assert_close(DD, 2.124343822e-07, 1e-5, '2pcf DD')
    assert_close(RR, 3.362627962e-07, 1e-5, '2pcf RR')


    out3 = 'tmp_3pcf.out'
    cmd3 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 6 -wgt -nmu 1 -out {out3} -3pcf"
    run(cmd3)
    tmp=np.loadtxt(out3,skiprows=1)

    DDD=np.mean(tmp[:,6])
    RRR=np.mean(tmp[:,7])

    assert_close(DDD, 5.960412001428572e-12, 1e-5, '3pcf DDD')
    assert_close(RRR, 5.04928030545e-12, 1e-5, '3pcf DDD')

    os.remove(out2)
    os.remove(out3)
    print('Correlation tests passed')

if __name__ == '__main__':
    main()
