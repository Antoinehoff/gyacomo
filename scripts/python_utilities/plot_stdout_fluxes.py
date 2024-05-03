import sys
import re
import matplotlib.pyplot as plt
import load_data as loader

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py filename")
        sys.exit(1)

    filename = sys.argv[1]
    dict = loader.read_data_std(filename)
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