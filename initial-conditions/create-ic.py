import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

def create_heatmap(name, cmap='inferno'):
    array = np.load(name+'.npy')
    plt.figure(figsize=(8, 6))
    plt.imshow(array, cmap=cmap, origin="lower",vmin=-1,vmax=1)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(name)
    plt.savefig(f'{name}.pdf', bbox_inches='tight')
    print(f'View new initial condition at {name}.pdf')

def create_uniform_mix(path, order_param, delta, dim=250, seed=12041997):
    
    '''
    path: where to save npy file 
    concentration: float in [-1,1]. -1 = completely salt, 1 = completely solvent
    delta: amplitude of noise to be added to concentration. 
    dim: side length of initial condition
    '''

    np.random.seed(seed)

    A = order_param*np.ones((dim, dim))
    noise = np.random.uniform(-delta, delta, size = (dim, dim))
    A+=noise

    np.save(path, A)
    create_heatmap(path)
    return A


if __name__ == "__main__":
    order_params = [0.2, 0.6] # Dynamics of Cahn Hilliard Eq are symmetric about order parameter = 0. No need to study negative order params. 0.2 -> fig_3, 0.6 -> fig_s22.
    deltas = [0.01]
    for order_param in order_params:
        for delta in deltas:
            path = str(SCRIPT_DIR / f'{order_param}-{delta}')
            create_uniform_mix(path=path, order_param=order_param, delta=delta)