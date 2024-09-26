import numpy as np
import matplotlib.pyplot as plt
import argparse

def parse_data(filename, selected_steps):
    """
    Parses the data from the input file and returns numpy arrays for x, E, phi, n, and nion for selected time steps.
    Assumes the input data has the following columns:
    x | nion | phi | n | E
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract data after the headers and ignore any 'Time Step' lines
    time_steps = []
    data_dict = {}
    current_step = None

    for line in lines:
        if "Time Step:" in line:
            # Extract time step number
            try:
                step_number = int(line.split(':')[1].strip())
                current_step = step_number
                data_dict[current_step] = {'x': [], 'nion': [], 'phi': [], 'n': [], 'E': []}
            except ValueError:
                current_step = None
        elif "------" in line:
            continue
        elif current_step is not None:
            try:
                values = line.split()
                if len(values) >= 5:
                    x = float(values[0])
                    E = float(values[1].replace('E+', 'e').replace('E-', 'e'))
                    phi = float(values[2])
                    n = float(values[3].replace('E+', 'e').replace('E-', 'e'))
                    nion = 0 if '*' in values[4] else float(values[4].replace('E+', 'e').replace('E-', 'e'))
                    
                    data_dict[current_step]['x'].append(x)
                    data_dict[current_step]['E'].append(E)
                    data_dict[current_step]['phi'].append(phi)
                    data_dict[current_step]['n'].append(n)
                    data_dict[current_step]['nion'].append(nion)
            except ValueError:
                # Skip lines with invalid data
                pass
    
    # Extract data for selected time steps
    x_all, E_all, phi_all, n_all, nion_all = [], [], [], [], []
    for step in selected_steps:
        if step in data_dict:
            x_all.append(np.array(data_dict[step]['x']))
            E_all.append(np.array(data_dict[step]['E']))
            phi_all.append(np.array(data_dict[step]['phi']))
            n_all.append(np.array(data_dict[step]['n']))
            nion_all.append(np.array(data_dict[step]['nion']))

    return x_all, E_all, phi_all, n_all, nion_all

def plot_data(x_all, E_all, phi_all, n_all, nion_all, selected_steps, colors):
    """
    Plots E, phi, n, and nion as functions of x for multiple time steps with specified colors.
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    for step, x, E, phi, n, nion, color in zip(selected_steps, x_all, E_all, phi_all, n_all, nion_all, colors):
        # Plot E vs x
        axes[0, 0].plot(x, E, label=f'Time Step: {step}', color=color, marker='o')
        axes[0, 0].set_xlabel('x (m)')
        axes[0, 0].set_ylabel('nion (1/m^3)')
        axes[0, 0].set_title('Ion Density (nion) vs x')
        axes[0, 0].grid(True)

        # Plot phi vs x
        axes[0, 1].plot(x, phi, label=f'Time Step: {step}', color=color, marker='o')
        axes[0, 1].set_xlabel('x (m)')
        axes[0, 1].set_ylabel('phi (V)')
        axes[0, 1].set_title('Electric Potential (phi) vs x')
        axes[0, 1].grid(True)

        # Plot n vs x
        axes[1, 0].plot(x, n, label=f'Time Step: {step}', color=color, marker='o')
        axes[1, 0].set_xlabel('x (m)')
        axes[1, 0].set_ylabel('n (1/m^3)')
        axes[1, 0].set_title('Electron Density (n) vs x')
        axes[1, 0].grid(True)

        # Plot nion vs x
        axes[1, 1].plot(x, nion, label=f'Time Step: {step}', color=color, marker='o')
        axes[1, 1].set_xlabel('x (m)')
        axes[1, 1].set_ylabel('E (V/m)')
        axes[1, 1].set_title('Electric Field (E) vs x')
        axes[1, 1].grid(True)

    # Add legends to all plots
    for ax in axes.flat:
        ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot data from simulation output file.')
    parser.add_argument('filename', type=str, help='Path to the data file')
    parser.add_argument('steps', type=int, nargs='+', help='List of time steps to plot')
    parser.add_argument('--colors', type=str, nargs='+', default=None, help='List of colors for each time step')

    args = parser.parse_args()

    # Ensure the number of colors matches the number of time steps
    if args.colors and len(args.colors) != len(args.steps):
        parser.error('The number of colors must match the number of time steps.')

    x_all, E_all, phi_all, n_all, nion_all = parse_data(args.filename, args.steps)
    
    # Use default colors if none are provided
    colors = args.colors if args.colors else plt.get_cmap('tab10').colors[:len(args.steps)]

    plot_data(x_all, E_all, phi_all, n_all, nion_all, args.steps, colors)
