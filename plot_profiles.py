import matplotlib.pyplot as plt
import argparse
import pandas as pd

# Set up the argument parser
p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument("-output_files", type=str, nargs='+', help="Space-separated list of data files to be plotted")
p.add_argument("-colors", type=str, nargs='+', help="Space-separated list of colors for the different data files")
args = p.parse_args()

# Read data from each file into a list of DataFrames
data = [pd.read_csv(f.strip(), sep='\s+', header=None, skiprows=2, 
                    names=["x", "nion", "phi", "n", "E"]) for f in args.output_files]

# Create a 2x2 subplot
fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

# Loop through each dataset and plot
for i, df in enumerate(data):
    color = args.colors[i % len(args.colors)]  # Cycle through colors if more files than colors
    filename = args.output_files[i].split('/')[-1]  # Get only the filename without the path
    
    # Plot nion vs x
    axes[0, 0].plot(df["x"], df["nion"], color=color, marker='o', label=filename)
    axes[0, 0].set_xlabel('x (m)')
    axes[0, 0].set_ylabel('nion (1/m^3)')
    axes[0, 0].set_title('Ion Density (nion) vs x')
    axes[0, 0].grid(True)
    axes[0, 0].legend()  # Add legend to this subplot

    # Plot phi vs x
    axes[0, 1].plot(df["x"], df["phi"], color=color, marker='o', label=filename)
    axes[0, 1].set_xlabel('x (m)')
    axes[0, 1].set_ylabel('phi (V)')
    axes[0, 1].set_title('Electric Potential (phi) vs x')
    axes[0, 1].grid(True)
    axes[0, 1].legend()  # Add legend to this subplot

    # Plot n vs x
    axes[1, 0].plot(df["x"], df["n"], color=color, marker='o', label=filename)
    axes[1, 0].set_xlabel('x (m)')
    axes[1, 0].set_ylabel('n (1/m^3)')
    axes[1, 0].set_title('Electron Density (n) vs x')
    axes[1, 0].grid(True)
    axes[1, 0].legend()  # Add legend to this subplot

    # Plot E vs x
    axes[1, 1].plot(df["x"], df["E"], color=color, marker='o', label=filename)
    axes[1, 1].set_xlabel('x (m)')
    axes[1, 1].set_ylabel('E (V/m)')
    axes[1, 1].set_title('Electric Field (E) vs x')
    axes[1, 1].grid(True)
    axes[1, 1].legend()  # Add legend to this subplot

# Display the plots
plt.show()