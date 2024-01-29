import sys
import numpy as np

def csv_maxdiff(csv1, csv2):
    reference = np.array(list(np.loadtxt(open(csv1), delimiter=","))).astype("float")
    subject = np.array(list(np.loadtxt(open(csv2), delimiter=","))).astype("float")
    diff = np.abs((reference-subject))
    maxdiff_indices = np.argmax(diff)
    maxdiff = diff[np.unravel_index(maxdiff_indices, diff.shape)]
    v1 = reference[np.unravel_index(maxdiff_indices, reference.shape)]
    v2 = subject[np.unravel_index(maxdiff_indices, subject.shape)]
    return maxdiff, v1, v2

def main():
    csv1 = sys.argv[1]
    csv2 = sys.argv[2]
    pass_limit = float(sys.argv[3])

    maxdiff, v1, v2 = csv_maxdiff(csv1, csv2)
    print(f'{csv1} and {csv2}: maxdiff {maxdiff:.8f} between {v1:.8f} and {v2:.8f}')

    if maxdiff > pass_limit:
        print('Fail')
        sys.exit(1)
    else:
        print('Pass')

if __name__ == "__main__":
    main()
