from pathlib import Path

from pyd4 import D4File

def main():
    test_file = Path("../../tests/testthat/testdata/test.d4")
    file = D4File(str(test_file))

    chroms = dict(file.chroms())

    # print(file[("chr1", 0, 1000)])
    # print(file.resample("chr1:0-1000", method="mean", bin_size=10))
    hist = file.histogram("chr1:0-1000", min=99, max=200)
    print(hist.std())




if __name__ == '__main__':
    main()