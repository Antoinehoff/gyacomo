import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

def load_transport_data(file_path):
    with h5py.File(file_path, 'r') as file:
        # Check if the group exists
        group_path = "/data/var0d"
        if group_path not in file:
            print(f"Group '{group_path}' not found in the HDF5 file.")
            return None

        group = file[group_path]
        array_data = {key: np.array(group[key]) for key in group.keys()}

    return array_data

def plot_time_traces(time, data, key):
    plt.plot(time, data, label=key)
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.show()

def read_namelist_from_h5(file_path, namelist_path):
    with h5py.File(file_path, 'r') as h5_file:
        # Check if the namelist path exists
        if namelist_path not in h5_file:
            print(f"Namelist path '{namelist_path}' not found in the HDF5 file.")
            return None

        namelist_group = h5_file[namelist_path]

        # Create a dictionary to store the namelist parameters
        namelist_params = {}

        # Iterate through items in the group
        for key, value in namelist_group.items():
            # Check if the item is a dataset (parameter)
            if isinstance(value, h5py.Dataset):
                namelist_params[key] = value[()]  # Add the parameter to the dictionary

    return namelist_params

# Check for the correct number of command-line arguments
if len(sys.argv) != 2:
    print("Usage: python script.py <file_path>")
    sys.exit(1)

# Get file path from command-line arguments
file_path = sys.argv[1]

namelist_params = read_namelist_from_h5(file_path, "/files/STDIN.00")

# Print the namelist parameters
if namelist_params is not None:
    print("Namelist Parameters:")
    for key, value in namelist_params.items():
        print(f"{key}: {value}")

# Load all data in the group /data/var0d
loaded_data = load_transport_data(file_path)
time        = loaded_data['time']
del loaded_data['time']
del loaded_data['cstep']

# Plot it
if loaded_data is not None:
    # Print all arrays in the loaded group
    print(f"Arrays loaded from '{file_path}' in group '/data/var0d':")
    for key, value in loaded_data.items():
        half_length = len(value) // 2
        trace = value[:half_length]
        plot_time_traces(time,trace,key)
