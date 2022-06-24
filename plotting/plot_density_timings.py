import numpy as np
import matplotlib.pyplot as plt

def main():
    data_fixed = np.loadtxt("density_timings_fixedh.csv", delimiter=',')
    N_fixed = data_fixed[:, 0]
    time_fixed = data_fixed[:, 1] * 1000

    data_var = np.loadtxt("density_timings_varh.csv", delimiter=',')
    N_var = data_var[:, 0]
    time_var = data_var[:, 1] * 1000

    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = 'lightgrey'
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'cm'
    plt.rcParams['font.size'] = 12

    fig = plt.figure(figsize=(2 * 3.7, 3))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(N_fixed, time_fixed, c='k')
    ax.set_xlabel("N")
    ax.set_ylabel("Time (ms)")
    ax.text(0.05, 0.9, "a)", transform=ax.transAxes)
    ax = fig.add_subplot(1, 2, 2)
    ax.plot(N_var, time_var, c='k')
    ax.set_xlabel("N")
    ax.text(0.05, 0.9, "b)", transform=ax.transAxes)
    plt.tight_layout()
    plt.savefig("density_timings.pdf")

if __name__ == "__main__":
    main()
