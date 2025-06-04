import os
import subprocess
import math

def run(cmd):
    print(cmd)
    subprocess.check_call(cmd, shell=True)

def average_from_file(path, col):
    values = []
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.split()
            if len(parts) > col:
                try:
                    val = float(parts[col])
                except ValueError:
                    continue
                if math.isfinite(val):
                    values.append(val)
    return sum(values) / len(values)

def assert_close(val, expected, tol, msg):
    if abs(val - expected) > tol:
        raise AssertionError(f"{msg}: got {val} expected {expected}")

def main():
    bindir = os.path.join('..', 'bin')
    gal = os.path.join('..', 'example', 'test.gal')
    ran = os.path.join('..', 'example', 'test.ran')

    out2 = 'tmp_2pcf.out'
    cmd2 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 10 -wgt -nmu 2 -RSD -out {out2} -2pcf"
    run(cmd2)
    mean2 = average_from_file(out2, 6)
    assert_close(mean2, 2.25823, 1e-5, '2pcf mean')

    out3 = 'tmp_3pcf.out'
    cmd3 = f"{os.path.join(bindir, 'gramsci')} -gal {gal} -ran {ran} -rmin 1.0 -rmax 30.0 -nbins 6 -wgt -nmu 2 -RSD -out {out3} -3pcf"
    run(cmd3)
    mean3 = average_from_file(out3, 10)
    assert_close(mean3, 8.18247, 1e-4, '3pcf mean')

    os.remove(out2)
    os.remove(out3)
    print('Correlation tests passed')

if __name__ == '__main__':
    main()
