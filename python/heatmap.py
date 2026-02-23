#!/usr/bin/env python

# general imports
import argparse as ap, matplotlib.animation as anm, matplotlib.pyplot as plt, matplotlib.ticker as tck, numpy as np

# additional imports
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# set rcParams
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Libertinus Serif"],
    "text.latex.preamble": r"\usepackage{braket}"
})

def heatmap():

    # create the parser
    parser = ap.ArgumentParser(
        prog="Zinq Heatmap Plotter", description="Heatmap plotting script for the Zinq package.",
        formatter_class=lambda prog: ap.HelpFormatter(prog, max_help_position=128),
        add_help=False, allow_abbrev=False
    )

    # add the arguments for the parser
    parser.add_argument("-h", "--help", action="help", default=ap.SUPPRESS, help="This help message.")
    parser.add_argument("-o", "--output", type=str, help="The output file to save the plot.")
    parser.add_argument("-a", "--animate", type=int, help="Animate the plot.")

    # add the plotting arguments
    parser.add_argument("--colormap", type=str, default="viridis", help="The colormap to use for the plot.")
    parser.add_argument("--dpi", type=int, default=96, help="The DPI of the plot.")
    parser.add_argument("--figsize", type=float, nargs=2, default=[4, 4], help="The dimensions of the plot in inches.")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second for the animation.")
    parser.add_argument("--title", type=str, nargs="+", help="The title of the plot.")
    parser.add_argument("--transform", type=str, default="sum", help="The transform function used when multiple columns are plotted.")
    parser.add_argument("--xlabel", type=str, nargs="+", help="The an x-axis label.")
    parser.add_argument("--xtickon", type=int, nargs="+", help="Whether to show the x-axis ticks for each subplot.")
    parser.add_argument("--xticklabelon", type=int, nargs="+", help="Whether to show the x-axis tick labels for each subplot.")
    parser.add_argument("--xlim", type=float, nargs="+", help="The x-axis limits of the plot.")
    parser.add_argument("--ylabel", type=str, nargs="+", help="The an y-axis label.")
    parser.add_argument("--ylim", type=float, nargs="+", help="The y-axis limits of the plot.")
    parser.add_argument("--ytickon", type=int, nargs="+", help="Whether to show the y-axis ticks for each subplot.")
    parser.add_argument("--yticklabelon", type=int, nargs="+", help="Whether to show the y-axis tick labels for each subplot.")

    # add positional arguments and parse the arguments
    parser.add_argument("files", nargs="+", help="The data files to plot."); args = parser.parse_args()

    # function to load the data from the input string as well as the columns to plot
    def load(fname): tmp = fname.split(":"); data = np.loadtxt(tmp[0], ndmin=2, skiprows=1); return data, np.array(tmp[1].split(","), int) if len(tmp) == 2 and tmp[1] else np.arange(data.shape[1] - 2)

    # function to return the iterator for the data and columns
    dcit = lambda d, c: ((d_i, j) for i, (d_i, c_i) in enumerate(zip(d, c)) for j in c_i)

    # load the data with plotted columns and get the total number of lines
    data, data_cols = zip(*[load(file) for file in args.files]); nline = sum(len(data_cols[i]) for i in range(len(data_cols)))

    # define the transform function
    if   args.transform == "sum":  transform = lambda z: np.sum        (z, axis=1)
    elif args.transform == "norm": transform = lambda z: np.linalg.norm(z, axis=1)
    else: raise ValueError(f"TRANSFORM '{args.transform}' NOT RECOGNIZED"        )

    # calculate the number of rows and columns for the plot
    nrow = [i for i in range(1, int(np.sqrt(len(args.files))) + 1) if len(args.files) % i == 0][-1]; ncol = len(args.files) // nrow

    # create the figure and the container for plots
    fig, ax = plt.subplots(nrow, ncol, dpi=args.dpi, figsize=(ncol * args.figsize[1], nrow * args.figsize[0])); plots = []

    # reshape the axes if the axes are not a matrix
    ax = np.array([ax]) if not isinstance(ax, np.ndarray) else (ax if not isinstance(ax[0], np.ndarray) else ax.flatten())

    # loop over data and its columns to plot
    for data_i, data_cols_i in zip(data, data_cols):

        # define the independent coordinates
        x = data_i[:, 0].reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))
        y = data_i[:, 1].reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))

        # define the dependent coordinates
        z = transform(data_i[:, data_cols_i + 2]).reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0])))

        # calculate the z limits
        zmin = np.min([transform(data_i[:, data_cols_i + i * (args.animate if args.animate else 0) + 2]).min() for i in range((data[0].shape[1] - 2) // args.animate if args.animate else 1)])
        zmax = np.max([transform(data_i[:, data_cols_i + i * (args.animate if args.animate else 0) + 2]).max() for i in range((data[0].shape[1] - 2) // args.animate if args.animate else 1)])

        # plot the data and append the plot to the list
        plots.append(ax[len(plots)].pcolormesh(x, y, z, cmap=args.colormap, edgecolor="face", linewidth=1e-4, vmin=zmin, vmax=zmax))

    # set the title
    if args.title: [ax[i].set_title(args.title[i]) for i in range(len(args.title)) if args.title[i]]

    # set the axis labels
    if args.xlabel: [ax[i].set_xlabel(args.xlabel[i]) for i in range(len(args.xlabel)) if args.xlabel[i]]
    if args.ylabel: [ax[i].set_ylabel(args.ylabel[i]) for i in range(len(args.ylabel)) if args.ylabel[i]]

    # set the axis limits
    if args.xlim: [ax[i].set_xlim([args.xlim[2 * i], args.xlim[2 * i + 1]]) for i in range(0, len(args.xlim) // 2) if not np.isnan(args.xlim[2 * i])]
    if args.ylim: [ax[i].set_ylim([args.ylim[2 * i], args.ylim[2 * i + 1]]) for i in range(0, len(args.ylim) // 2) if not np.isnan(args.ylim[2 * i])]

    # loop over the axes
    for i in range(len(ax)):

        # enable or disable the x-axis and y-axis ticks
        if args.xtickon and len(args.xtickon) > i: ax[i].xaxis.set_tick_params(bottom=args.xtickon[i])
        if args.ytickon and len(args.ytickon) > i: ax[i].yaxis.set_tick_params(left  =args.ytickon[i])

        # enable or disable the x-axis tick labels
        if args.xticklabelon and len(args.xticklabelon) > i: ax[i].xaxis.set_tick_params(labelbottom=args.xticklabelon[i])
        if args.yticklabelon and len(args.yticklabelon) > i: ax[i].yaxis.set_tick_params(labelleft  =args.yticklabelon[i])

    # set the layout
    fig.tight_layout()

    # update function
    def update(frame):

        # loop over data and its columns to plot
        for i, (data_i, data_cols_i) in enumerate(zip(data, data_cols)):

            # update the data
            plots[i].set_array(transform(data_i[:, data_cols_i + frame * args.animate + 2]).reshape(int(np.sqrt(data_i.shape[0])), int(np.sqrt(data_i.shape[0]))))

    # create the animation
    anim = anm.FuncAnimation(fig, update, frames=(data[0].shape[1] - 2) // args.animate, init_func=lambda: None, interval=1000 // args.fps) if args.animate else None

    # set the window title
    if not args.output: fig.canvas.manager.set_window_title("Zinq Heatmap Plotter")

    # save or show the plot
    (anim.save(args.output, fps=args.fps) if args.animate else plt.savefig(args.output)) if args.output else plt.show()

if __name__ == "__main__": heatmap()
