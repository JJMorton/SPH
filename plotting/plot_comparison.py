#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from os import path

def main(datafiles):

    print("Reading data...")

    # Read the metadata from the file header of the first file
    metadata = {}
    with open(datafiles[0]) as f:
        firstline = next(f).replace('\n', '')
        kv = [ tuple(var.split("=")) for var in firstline.split("\t") ]
        metadata_str = { var[0]: var[1] for var in kv }
        for k, v in metadata_str.items():
            try:
                metadata[k] = np.round(float(v), 4)
            except ValueError:
                metadata[k] = v
    print(metadata)

    data = np.loadtxt(datafiles[0], skiprows=1)
    position_first = data[0::6, :]
    density_first = data[1::6, :]
    velocity_first = data[2::6, :]
    internal_first = data[3::6, :]
    pressure_first = data[5::6, :]
    print(position_first.shape)

    position_other = []
    density_other = []
    internal_other = []
    for datafile in datafiles[1:]:
        data = np.loadtxt(datafile, skiprows=1)
        position_other.append(data[0::6, :])
        density_other.append(data[1::6, :])
        internal_other.append(data[3::6, :])
        print(position_other[-1].shape)

    # Non-dimensionalise
    if metadata["init"] == "colliding_streams":
        scale = density_first[0, 0]
        density_first /= scale
        pressure_first /= scale
        for density in density_other: density /= density[0, 0]
    if metadata["init"] == "shock_tube":
        print("[shock_tube] Density LHS: {:.3f}, density RHS: {:.3f}, ratio={:.3f}".format(density_first[0, 0], density_first[0, -1], density_first[0, 0] / density_first[0, -1]))
        gamma = metadata["gamma"]
        # internal *= gamma
        # velocity *= np.sqrt(gamma)
        scale = density_first[0, 0]
        density_first /= scale
        pressure_first /= scale
        for density in density_other: density /= density[0, 0]

    time = np.linspace(0, metadata["dump_interval"] * position_first.shape[0], position_first.shape[0])

    def calc_limits(arr, margin=0.1):
        lim = [np.amin(arr), np.amax(arr)]
        lim[0] -= (lim[1] - lim[0]) * margin
        lim[1] += (lim[1] - lim[0]) * margin
        return lim

    position_lim = (-metadata["box_width"]/2, metadata["box_width"]/2)

    # Calculate analytical solutions
    x_analytic = density_analytic = internal_analytic = pressure_analytic = velocity_analytic = None;
    if metadata["init"] == "shock_tube" and metadata["EOS"] == "ideal":
        data = np.loadtxt(path.join(sys.path[0], "shock_tube_analytical.csv"), delimiter=',', skiprows=1)
        x_analytic = data[:, 0]
        density_analytic = data[:, 1]
        internal_analytic = data[:, 2]
        pressure_analytic = data[:, 3]
        velocity_analytic = data[:, 4]

    elif metadata["init"] == "colliding_streams":
        timeend = time[-1]
        v0 = np.amax(velocity_first[0])
        rho0 = density_first[0, 0]
        p0 = pressure_first[0, 0]
        u0 = internal_first[0, 0]
        x_analytic = np.linspace(-metadata["box_width"]/2, metadata["box_width"]/2, 1000)
        x_centre = x_analytic[len(x_analytic) // 2]
        cs = v = rho1 = p1 = u1 = 0.0
        if metadata["EOS"] == "ideal":
            gamma = metadata["gamma"]
            cs = np.sqrt(gamma)
            v = 0.25 * ((gamma - 3)*v0 + np.sqrt((gamma + 1)**2 * v0**2 + 16*gamma))
            rho1 = rho0 * (1 + v0/v)
            p1 = p0 * (1 + v0 * (v + v0))
            u1 = u0 * (rho0  / rho1 + v0*v)
        elif metadata["EOS"] == "isothermal":
            cs = 1
            v = 0.5 * (-v0 + np.sqrt(v0**2 + 4))
            rho1 = rho0 * (1 + v0/v)
            p1 = p0 * (1 + v0/v)
            u1 = u0
        print("Calculated colliding streams analytical form with: v0 = {:.3f}, v = {:.3f}, cs = {:.3f}".format(v0, v, cs))
        position_lim = (-v * timeend * 2, v * timeend * 2)
        density_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, rho1, rho0)
        pressure_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, p1, p0)
        internal_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, u1, u0)
        velocity_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, 0, v0)
        velocity_analytic[len(velocity_analytic)//2:] *= -1

    print("Plotting...")
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = 'cm'
    plt.rcParams['font.size'] = 12
    num_plots = 1 + int(metadata["EOS"] == "ideal")
    plot_layout = (1, num_plots) # (rows, cols)
    fig = plt.figure(figsize=(3 * plot_layout[1], 3 * plot_layout[0]))
    gs = plt.GridSpec(*plot_layout, figure=fig)
    axes = list(gs.subplots().flatten())

    def plot_box(ax, box_width):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

    def plot_property(ax, pos, prop, ls):
        # return ax.plot(pos[0], prop[0], c="black", lw=1, ls='', marker='.', ms=0.5)
        ret = ax.plot(pos[-1], prop[-1], c="black", lw=1, ls=ls)
        ax.set_xlabel("Position")
        ax.set_xlim(*position_lim)
        ax.set_ylim(*calc_limits(prop))
        plot_box(ax, metadata["box_width"])
        return ret

    linestyles = ('-.', '--', ':')
    linestyle = 0

    # Plot density dist
    ax_density = axes.pop(0)
    if density_analytic is not None:
        ax_density.plot(x_analytic, density_analytic, c="lightblue", lw=4)
    plot_property(ax_density, position_first, density_first, ls='-')
    for position, density in zip(position_other, density_other):
        plot_property(ax_density, position, density, ls=linestyles[linestyle])
        linestyle += 1
    ax_density.set_ylabel("Density")

    linestyle = 0

    # Plot internal energy
    if metadata["EOS"] == "ideal":
        ax_internal = axes.pop(0)
        if internal_analytic is not None:
            ax_internal.plot(x_analytic, internal_analytic, c="lightblue", lw=4)
        plot_property(ax_internal, position_first, internal_first, ls='-')
        for position, internal in zip(position_other, internal_other):
            plot_property(ax_internal, position, internal, ls=linestyles[linestyle])
            linestyle += 1
        ax_internal.set_ylabel("Internal Energy")

    for ax in axes: ax.set_visible(False)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Provide data file")
        exit(1)

    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = 'lightgrey'
    plt.rcParams['font.size'] = 11

    main(sys.argv[1:])

