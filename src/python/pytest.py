from pathlib import Path

from pyd4 import D4File

def main():
    test_file = Path("../../tests/testthat/testdata/test.d4")
    file = D4File(str(test_file))

    chroms = dict(file.chroms())

    # print(file[("chr1", 0, 1000)])
    # print(file.resample("chr1:0-1000", method="mean", bin_size=10))
    hist = file.histogram("chr1:0-1000", min=99, max=200)

    print(f"median depth of region: {file.median(('chr1', 12, 22))}")
    print(f"below {hist.below}")
    print(f"above: {hist.above}")
    print(f"mean: {hist.mean()}")
    print(f"median: {hist.median()}")
    print(f"percentile: {hist.percentile(99.79)}")
    print(f"fraction_below: {hist.percentile_below(199)}")
    print(f"value_count(100): {hist.value_count(100)}")
    print(f"value_percentile(100): {hist.value_percentage(100)}")
    print(f"std: {hist.std()}")




if __name__ == '__main__':
    main()