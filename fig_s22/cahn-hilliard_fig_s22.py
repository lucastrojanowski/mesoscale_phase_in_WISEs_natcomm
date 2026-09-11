import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.colorbar as mcolorbar
import imageio as iio
import time
import sys
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['font.size']=14

def progressbar(it, size=60, out=sys.stdout):
    print('Working ... ')
    count = len(it)
    start = time.time()
    def show(j):
        x = int(size*j/count)

        remaining = ((time.time() - start) / j) * (count - j)        
        mins, sec = divmod(remaining, 60) # limit to minutes
        time_str = f"{int(mins):02}:{sec:03.1f}"
        print(f"Est wait {time_str}", end='\r', file=out, flush=True)
    show(0.1) # avoid div/0 
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)

def solve_cahn_hilliard(diff_coeff, dt, N, grid_size, save_path, gamma, snapshot_steps, u0=None, mp4_freq=100):
    """
    Solves the 2D Cahn-Hilliard equation with periodic boundary conditions and generates a heatmap MP4 video using imageio.
    Also extracts the frames listed in `snapshot_steps` and saves each as its own labeled-free snapshot figure (PDF and PNG).

    Parameters:
    - diff_coeff: Diffusion coefficient (float)
    - dt: Time step (float)
    - N: Number of time steps (int)
    - grid_size: Size of the 2D grid (int, e.g., 100 for a 100x100 grid)
    - save_path: Path to save the generated MP4 file (str)
    - gamma: Parameter in the Cahn-Hilliard equation (float)
    - snapshot_steps: Time step indices (t) to extract for the individual snapshot figures
    - u0: Initial condition (2D numpy array)
    """

    # Initialize parameters
    dx = 1  # Spatial resolution in x direction
    dy = 1  # Spatial resolution in y direction
    nx, ny = grid_size, grid_size  # Number of grid points in x and y
    nt = N  # Number of time steps
    alpha_x = diff_coeff*dt/(dx**2)  # Stability condition in x direction
    alpha_y = diff_coeff*dt/(dy**2)  # Stability condition in y direction

    # Check stability
    if (gamma*alpha_x > 0.01) or (gamma*alpha_y > 0.01):
        raise ValueError(f'Dγdt/(dx^2) = {gamma*alpha_x} !<= 0.01. Change parameters for stable simulations.')
    if (alpha_x > 0.1) or (alpha_y > 0.1):
        raise(ValueError(f'Ddt/(dx^2) = {alpha_x}!<= 0.1. Change parameters for stable simulations.'))

    if u0 is None:
        u0 = 0.5*np.ones((nx, ny))+0.1*np.random.uniform(0, 1, size=(nx, ny))

    u = u0.copy()

    # Compute Laplacian
    def laplacian(field):
        return (np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0) +
                np.roll(field, 1, axis=1) + np.roll(field, -1, axis=1) -4*field)/(dx*dy)

    # MP4 Writer
    writer = iio.get_writer(save_path, fps=100, codec="libx264", quality=10)
    fig_width = 1280
    fig_height = 720
    dpi = 80
    fig_size = (fig_width / dpi, fig_height / dpi)

    snapshots = {}  # t -> (volume fraction field, time value)

    for t in progressbar(range(nt+1)):
        lap_u = laplacian(u)
        lap_u3 = laplacian(u**3-u) #Bulk energy term
        lap_lap_u = laplacian(lap_u) # Interfacial energy term
        u+=dt*diff_coeff*(lap_u3-gamma*lap_lap_u)

        if t%mp4_freq == 0:
            volume_fraction = 1/2*(u+1) # u has asymptotic values between -1 and 1. Shift this to 0 and 1 by plotting 1/2(u+1). Now this can be interpreted as a volume fraction.

            fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
            ax.set_title(f"t={t*dt:.2f}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            heatmap = ax.imshow(volume_fraction, cmap="plasma", origin="lower", vmin=0, vmax=1)

            plt.colorbar(heatmap, label="Local Solvent Volume Fraction")
            fig.canvas.draw()
            frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
            frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))[..., :3]
            writer.append_data(frame)
            plt.close(fig)

            if t in snapshot_steps:
                snapshots[t] = (volume_fraction.copy(), t*dt)

    writer.close()
    print(f"MP4 video saved at {save_path}")

    save_snapshot_figures(snapshots, snapshot_steps, save_path)


def save_snapshot_figures(snapshots, snapshot_steps, save_path):

    """
    Saves a PDF and PNG for each requested snapshot frame, keeping the axes box and
    tick marks but with no tick labels, axis labels, or title, plus a single standalone
    colorbar PDF and PNG (axis/ticks kept, no tick labels), and a side-by-side combined
    comparison PDF and PNG of all frames with fully labeled axes/colorbar and each panel
    labeled by its time, all alongside the MP4 (same directory/basename as save_path).
    """

    ordered_steps = [t for t in snapshot_steps if t in snapshots]
    if not ordered_steps:
        print("No snapshot frames were captured; skipping snapshot figures.")
        return

    base, _ = os.path.splitext(save_path)
    out_dir = os.path.dirname(base)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for t in ordered_steps:
        field, time_val = snapshots[t]
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(field, cmap="plasma", origin="lower", vmin=0, vmax=1)
        ax.tick_params(labelbottom=False, labelleft=False)
        frame_base = f"{base}_t{t}"
        fig.savefig(frame_base + '.pdf', bbox_inches='tight', pad_inches=0)
        fig.savefig(frame_base + '.png', bbox_inches='tight', pad_inches=0, dpi=300)
        plt.close(fig)
        print(f"Snapshot figure saved at {frame_base}.pdf and {frame_base}.png")

    # Standalone colorbar, axis/ticks kept but no tick labels
    fig_cb = plt.figure(figsize=(1.5, 6))
    ax_cb = fig_cb.add_axes([0.4, 0.05, 0.2, 0.9])
    norm = mcolors.Normalize(vmin=0, vmax=1)
    cb = mcolorbar.ColorbarBase(ax_cb, cmap=plt.get_cmap("plasma"), norm=norm)
    cb.set_ticks([0, 0.25, 0.5, 0.75, 1])
    cb.ax.tick_params(labelleft=False, labelright=False)
    cb_base = f"{base}_colorbar"
    fig_cb.savefig(cb_base + '.pdf', bbox_inches='tight', pad_inches=0)
    fig_cb.savefig(cb_base + '.png', bbox_inches='tight', pad_inches=0, dpi=300)
    plt.close(fig_cb)
    print(f"Colorbar figure saved at {cb_base}.pdf and {cb_base}.png")

    save_combined_figure(snapshots, ordered_steps, save_path)


def save_combined_figure(snapshots, ordered_steps, save_path):

    """
    Saves a single side-by-side comparison figure containing all snapshot frames,
    each panel directly labeled by its associated time, with x tick labels, y tick
    labels, the y axis, and a shared, fully labeled colorbar (with tick labels).
    """

    base, _ = os.path.splitext(save_path)
    n_panels = len(ordered_steps)

    fig, axes = plt.subplots(1, n_panels, figsize=(4.5*n_panels, 5), sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, t in zip(axes, ordered_steps):
        field, time_val = snapshots[t]
        heatmap = ax.imshow(field, cmap="plasma", origin="lower", vmin=0, vmax=1)
        ax.set_xlabel("x")
        ax.text(0.05, 0.95, f"t={time_val:.2f}", transform=ax.transAxes,
                ha='left', va='top', color='white', fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='black', alpha=0.5, pad=0.3))

    axes[0].set_ylabel("y")

    fig.colorbar(heatmap, ax=axes, label="Local Solvent Volume Fraction",
                 fraction=0.046, pad=0.04)

    combined_base = f"{base}_combined"
    fig.savefig(combined_base + '.pdf', bbox_inches='tight', pad_inches=0.1)
    fig.savefig(combined_base + '.png', bbox_inches='tight', pad_inches=0.1, dpi=300)
    plt.close(fig)
    print(f"Combined comparison figure saved at {combined_base}.pdf and {combined_base}.png")

# Example usage
if __name__ == '__main__':

    """
    Specify initial parameters:

    Parameters:
    
    - dt: Time step (float)
    - N: Number of time steps (int)
    - grid_size: Size of the 2D grid (int, e.g., 100 for a 100x100 grid)
    - diff_coeff: Diffusion coefficient (float)
    - save_path: Path to save the generated MP4 file (str)
    - gamma: Parameter in the Cahn-Hilliard equation (float)
    - u0: Initial condition (2D numpy array)
    - snapshot_frames: Time step indices (t) to extract for the individual snapshot figures

    """
    # Specifies simulation parameters

    dt = 0.01 # Integration step size
    N = 100001 # Number of steps
    grid_size = 250 # Square grid dimension - must match loaded initial condition numpy dimension
    diff_coeff = 1 # Diffusion coefficient
    gammas = [1] # Boundary energy penalty
    mp4_freq = 100 #Output every nth frame to mp4 file

    order_params = [0.6] # Take (order_param + 1)/2 to get solvent ratio. Fig 3: order_param = 0.2. Fig S22: order_param = 0.6

    # Time step indices (t) to extract as individual snapshot figures (PDF/PNG)

    snapshot_frames = (0, 7100, 100000)

    # Perform simulations of specified systems

    for gamma in gammas:
        for order_param in order_params:
            t1 = time.time()
            u0 = np.load(SCRIPT_DIR.parent / 'initial-conditions' / f'{order_param}-0.01.npy')
            save_path = str(SCRIPT_DIR / f'figs22_{round(1/2*(order_param+1),2)}.mp4')
            solve_cahn_hilliard(diff_coeff, dt, N, grid_size, save_path, gamma, snapshot_frames, u0=u0, mp4_freq = mp4_freq)
            t2 = time.time()
            print(f'Total execution time: {t2-t1} seconds')