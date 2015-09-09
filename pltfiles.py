import sys
import matplotlib.pyplot as plt

if (len(sys.argv) == 1):
    sys.exit("argv[1] = file name")
elif (len(sys.argv) == 2):
    x = []
    y = []

    fileName = sys.argv[1]
    try:
	ifile = open(fileName, 'r')
    except:
	sys.exit(fileName + " does not exist.")

    for index, string in enumerate(ifile):
	a = string.split()
	x.append(float(a[0]))
	y.append(float(a[1]))
    ifile.close()

    lowerY = 0
    if (min(y) < 0):
	lowerY = min(y)

    plt.plot(x, y, '-o')
    plt.xlim([0, max(x)])
    plt.ylim(min(y) - 0.2, 0.2+max(y))
    plt.grid()
    plt.show()

elif(len(sys.argv) == 3):
    fileOne = sys.argv[1]
    fileTwo = sys.argv[2]
    try:
	readFileOne = open(fileOne, 'r')
    except:
	sys.exit(fileOne + " does not exist. ")
    try:
	readFileTwo = open(fileTwo, 'r')
    except:
	sys.exit(fileTwo + " does not exist. ")
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    for index, string in enumerate(readFileOne):
	a = string.split()
	x1.append(float(a[0]))
	y1.append(float(a[1]))
    readFileOne.close()

    for index, string in enumerate(readFileTwo):
	a = string.split()
	x2.append(float(a[0]))
	y2.append(float(a[1]))
    readFileTwo.close()


    lowerY = 0
    if (min(y1) < 0 or min(y2) < 0):
	lowerY = min(min(y1), min(y2))

    plt.plot(x1, y1, 'b-o', x2, y2, 'r-o')
    plt.xlim([0, max(max(x1), max(x2))])
    plt.ylim([lowerY - 0.2, 0.2 + 1.2*max(max(y1), max(y2))])
    plt.grid()
    plt.show()
else:
    sys.exit("Only support one or two files. ")
