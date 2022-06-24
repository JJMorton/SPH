#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from os import path
from subprocess import run as run_process
from matplotlib.patches import Rectangle

def main(datafile):

    show_ghost = False
    makegif = False
    use_tex = False

    giffile = path.splitext(datafile)[0] + ".gif"
    pdffile = path.splitext(datafile)[0] + ".pdf"

    print("Reading data...")

    # Read the metadata from the file header
    metadata = {}
    with open(datafile) as f:
        firstline = next(f).replace('\n', '')
        kv = [ tuple(var.split("=")) for var in firstline.split("\t") ]
        metadata_str = { var[0]: var[1] for var in kv }
        for k, v in metadata_str.items():
            try:
                metadata[k] = np.round(float(v), 4)
            except ValueError:
                metadata[k] = v
    print(metadata)
    metadata_title = ""
    for i in range(len(metadata)):
        if i > 0:
            metadata_title += ", "
            if i % 5 == 0:
                metadata_title += '\n'
        metadata_title += "{}={}".format(*list(metadata.items())[i])

    data = np.genfromtxt(datafile, skip_header=1, filling_values=0, usecols=np.arange(int(metadata["N"])), delimiter='\t')
    # data = np.loadtxt(datafile, skiprows=1)
    position = data[0::6, :]
    density = data[1::6, :]
    velocity = data[2::6, :]
    internal = data[3::6, :]
    slength = data[4::6, :]
    pressure = data[5::6, :]

    # Non-dimensionalise
    if metadata["init"] == "colliding_streams":
        scale = density[0, 0]
        density /= scale
        pressure /= scale
    if metadata["init"] == "shock_tube":
        print("[shock_tube] Density LHS: {:.3f}, density RHS: {:.3f}, ratio={:.3f}".format(density[0, 0], density[0, -1], density[0, 0] / density[0, -1]))
        gamma = metadata["gamma"]
        # internal *= gamma
        # velocity *= np.sqrt(gamma)
        scale = density[0, 0]
        density /= scale
        pressure /= scale

    kinetic_total = 0.5 * np.sum(velocity * velocity, axis=1)
    internal_total = np.sum(internal, axis=1)
    momentum = np.sum(velocity, axis=1)
    time = np.linspace(0, metadata["dump_interval"] * position.shape[0], position.shape[0])

    def calc_limits(arr, margin=0.1):
        lim = [np.amin(arr), np.amax(arr)]
        lim[0] -= (lim[1] - lim[0]) * margin
        lim[1] += (lim[1] - lim[0]) * margin
        return lim

    position_lim = (-metadata["box_width"], metadata["box_width"]) if show_ghost else (-metadata["box_width"]/2, metadata["box_width"]/2)

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
        v0 = np.amax(velocity[0])
        rho0 = density[0, 0]
        p0 = pressure[0, 0]
        u0 = internal[0, 0]
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
        print("[colliding_streams] v0 = {:.3f}, v = {:.3f}, cs = {:.3f}".format(v0, v, cs))
        position_lim = (-v * timeend * 2, v * timeend * 2)
        density_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, rho1, rho0)
        pressure_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, p1, p0)
        internal_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, u1, u0)
        velocity_analytic = np.where(np.abs(x_analytic - x_centre) <= v * timeend, 0, v0)
        velocity_analytic[len(velocity_analytic)//2:] *= -1

    print("Plotting...")
    if use_tex:
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = 'cm'
        plt.rcParams['font.size'] = 12
    num_plots = 5 + int(metadata["EOS"] == "ideal") + int(metadata["var_smoothlength"])
    plot_layout = (3 if num_plots > 6 else 2, 3) # (rows, cols)
    fig = plt.figure(figsize=((4 if show_ghost else 3) * plot_layout[1], 3 * plot_layout[0]))
    gs = plt.GridSpec(*plot_layout, figure=fig)
    fig.suptitle(metadata_title, fontproperties={'family': 'monospace', 'size': 11}, c='grey')
    axes = list(gs.subplots().flatten())

    def plot_box(ax, box_width):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.add_patch(Rectangle((box_width/2, ylim[0]), xlim[1] - box_width/2, ylim[1] - ylim[0], facecolor='k', alpha=0.1, hatch='/', edgecolor='grey'))
        ax.add_patch(Rectangle((xlim[0], ylim[0]), -box_width/2 - xlim[0], ylim[1] - ylim[0], facecolor='k', alpha=0.1, hatch='/', edgecolor='grey'))

    def plot_property(ax, pos, prop):
        # return ax.plot(pos[0], prop[0], c="black", lw=1, ls='', marker='.', ms=0.5)
        ret = ax.plot(pos[-1], prop[-1], c="black", lw=1)
        ax.set_xlabel("Position")
        ax.set_xlim(*position_lim)
        ax.set_ylim(*calc_limits(prop))
        plot_box(ax, metadata["box_width"])
        return ret

    # Plot density dist
    ax_density = axes.pop(0)
    if density_analytic is not None:
        ax_density.plot(x_analytic, density_analytic, c="lightblue", lw=4)
    density_line, = plot_property(ax_density, position, density)
    ax_density.set_ylabel("Density")

    # Plot internal energy
    if metadata["EOS"] == "ideal":
        ax_internal = axes.pop(0)
        if internal_analytic is not None:
            ax_internal.plot(x_analytic, internal_analytic, c="lightblue", lw=4)
        internal_line, = plot_property(ax_internal, position, internal)
        ax_internal.set_ylabel("Internal Energy")

    # Plot pressure
    ax_pressure = axes.pop(0)
    if pressure_analytic is not None:
        ax_pressure.plot(x_analytic, pressure_analytic, c="lightblue", lw=4)
    pressure_line, = plot_property(ax_pressure, position, pressure)
    ax_pressure.set_ylabel("Pressure")

    # Plot position against velocity
    ax_phase = axes.pop(0)
    if velocity_analytic is not None:
        ax_phase.plot(x_analytic, velocity_analytic, c="lightblue", lw=4)
    phase_line, = plot_property(ax_phase, position, velocity)
    ax_phase.set_ylabel("Velocity")

    if metadata["var_smoothlength"]:
        # Plot smoothing length
        ax_slength = axes.pop(0)
        slength_line, = plot_property(ax_slength, position, slength)
        ax_slength.set_ylabel("Smoothing length")

    # Plot total energy as a function of time
    ax_energy = axes.pop(0)
    ax_energy.plot(time, kinetic_total, label="Kinetic", ls="--", c="black", lw=1)
    ax_energy.plot(time, internal_total, label="Internal", ls="-.", c="black", lw=1)
    ax_energy.plot(time, kinetic_total + internal_total, label="Sum", c="black", lw=1)
    ax_energy.set_ylabel("Total Energy")
    ax_energy.set_xlabel("Time")
    ax_energy.legend()

    # Plot total momentum as a function of time
    ax_mom = axes.pop(0)
    ax_mom.plot(time, momentum, c="black", lw=1)
    ax_mom.set_ylabel("Total momentum")
    ax_mom.set_xlabel("Time")
    momentum_lim = (min(np.amin(momentum), -1e-10), max(np.amax(momentum), 1e-10))
    ax_mom.set_ylim(*momentum_lim)

    for ax in axes: ax.set_visible(False)

    plt.tight_layout()
    fig.savefig(pdffile)
    if not makegif: return

    time_text = ax_density.text(0.0, 1.05, "t = {:.3f}".format(time[0]), transform=ax_density.transAxes, in_layout=True)
    plt.tight_layout()

    def animate(index):
        density_line.set_xdata(position[index])
        density_line.set_ydata(density[index])
        if metadata["EOS"] == "ideal":
            internal_line.set_xdata(position[index])
            internal_line.set_ydata(internal[index])
        pressure_line.set_xdata(position[index])
        pressure_line.set_ydata(pressure[index])
        phase_line.set_xdata(position[index])
        phase_line.set_ydata(velocity[index])
        if metadata["var_smoothlength"]:
            slength_line.set_xdata(position[index])
            slength_line.set_ydata(slength[index])
        time_text.set_text("[ t = {:.3f} ]".format(time[index]))

    frame_index = np.concatenate([np.array([0]*20), np.arange(position.shape[0]), np.array([-1]*20)])
    duration_per_unit_time = 5_000 # Number of animation ms per unit of simulation time
    animation = ani.FuncAnimation(fig, animate, frames=frame_index, interval=metadata["dump_interval"] * duration_per_unit_time, blit=True)

    print("Rendering animation...")
    animation.save(giffile)
    plt.close()
    print("Done")

    run_process(["xdg-open", giffile])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Provide data file")
        exit(1)

    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = 'lightgrey'
    plt.rcParams['font.size'] = 11

    main(sys.argv[1])

