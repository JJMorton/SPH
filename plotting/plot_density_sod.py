import numpy as np
import matplotlib.pyplot as plt

def main():
    data_var = np.loadtxt("density_sod_varh.csv", delimiter=',', skiprows=1)
    pos_var = data_var[:, 0]
    rho_var = data_var[:, 1]
    rho_var = rho_var / rho_var[0]

    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = 'lightgrey'
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'cm'
    plt.rcParams['font.size'] = 12

    fig = plt.figure(figsize=(3.7, 3))
    gs = plt.GridSpec(5, 1, figure=fig)
    ax = fig.add_subplot(gs[:4])
    ax.plot(pos_var, rho_var, c='k')
    ax.set_ylabel("Density")
    # ax.get_xaxis().set_visible(False)
    ax.set_xticklabels([])
    ax = fig.add_subplot(gs[4])
    ax.plot(pos_var, [0]*len(pos_var), c='k', ls='', marker='.', ms=1)
    ax.set_xlabel("Position")
    ax.set_yticks([])
    ax.grid(False)
    plt.tight_layout()
    # plt.subplots_adjust(hspace=0)
    plt.savefig("density_sod.pdf")

if __name__ == "__main__":
    main()
