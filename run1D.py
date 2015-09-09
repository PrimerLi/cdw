#!/usr/bin/env python 

def main():
    import os
    import sys

    T = []
    for i in range(100):
        T.append(0.02*i)
        os.system("./cdw1D.out " + str(T[i]))

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
