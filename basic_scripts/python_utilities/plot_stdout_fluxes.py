import sys
import re
import matplotlib.pyplot as plt

def read_data(filename):
    t_values = []
    Pxi_values = []
    Qxi_values = []
    dict       = {"t":[],"Pxi":[],"Pxe":[],"Qxi":[],"Qxe":[]}
    with open(filename, 'r') as file:
        for line in file:
            a = line.split('|')
            for i in a[1:-1]:
                b = i.split('=')
                dict[b[0].strip()].append(float(b[1]))

    return dict

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py filename")
        sys.exit(1)

    filename = sys.argv[1]
    dict = read_data(filename)
    time = dict["t"]
    dict.pop('t',None)

    plt.subplot(1, 2, 1)
    if dict["Pxi"]:
        plt.plot(time, dict["Pxi"],linestyle='-',color='red',label='ions')
    if dict["Pxe"]:
        plt.plot(time, dict["Pxe"],linestyle='-',color='blue',label='electrons')
    plt.xlabel('t')
    plt.ylabel('Px')
    plt.legend()

    plt.subplot(1, 2, 2)
    if dict["Qxi"]:
        plt.plot(time, dict["Qxi"],linestyle='-',color='red',label='ions')
    if dict["Qxe"]:
        plt.plot(time, dict["Qxe"],linestyle='-',color='blue',label='electrons')
    plt.xlabel('t')
    plt.ylabel('Qx')
    plt.legend()

    plt.tight_layout()
    plt.show()