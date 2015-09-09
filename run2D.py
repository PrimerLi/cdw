#!/usr/bin/env python 

def main():
    import os
    import sys

    T = []
    for i in range(150):
        T.append(0.01*i)
        os.system("./cdw2D.out " + str(T[i]))

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
