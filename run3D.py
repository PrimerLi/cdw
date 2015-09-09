#!/usr/bin/env python 

def main():
    import os
    import sys

    T = []
    for i in range(13):
        T.append(0.05*i)
        os.system("./cdw3D.out " + str(T[i]))

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
