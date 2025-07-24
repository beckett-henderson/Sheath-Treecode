#!/usr/bin/env python3
# plasma_sheath.py

"""
1D Plasma Sheath Simulation using Barnes–Hut Treecode
Author: Beckett Henderson
Date Created: 07/01/25
Date Last Modified: 07/21/25
This script simulates a 1D plasma sheath with immobile ions and mobile electrons.

It uses nondimensionalized variables based on the following reference scales:
    * Length scale: L = λ_D (Debye length)
    * Thermal speed: v_th (electron thermal velocity)
    * Time scale: T = L / v_th
    * Potential scale: φ_0 = m_e v_th^2 / q_0
With dimensionless variables:
    * x' = x / L,
    * v' = v / v_th,
    * t' = t / T,
    * φ' = φ / φ_0,
    * E' = (q_0 L) / (m_e v_th^2) E.
Time-stepping uses a leapfrog-like update:
    * v'_{n+1} = v'_n + E'(x') Δt',
    * x'_{n+1} = x'_n + v'_{n+1} Δt'.

All core parameters can be set via command-line flags; defaults match the author's reference setup.

Equations summary:
    * x = L x', v = v_th v', t = (L/v_th) t', φ = φ_0 φ'  (Scale definitions)
    * E' = (q_0 L)/(m_e v_th^2) E  (Field scaling)
    * dx'/dt' = v'  (Dimensionless kinematic equation)
    * dv'/dt' = E'  (Dimensionless dynamic equation)
    * ∂²φ'/∂x'² = -ρ'  (Dimensionless Poisson equation)
    * G(x',x'_p) = { x'(1–x'_p) if x'<x'_p; x'_p(1–x') otherwise }  (Dirichlet Green’s function)
    * A' = Σ q'_i G_left(x'_i), B' = Σ q'_i G_right(x'_i)  (Boundary correction terms)
    * v'_{n+1} = v'_n + E'(x') Δt', x'_{n+1} = x'_n + v'_{n+1} Δt'  (Leapfrog update)
"""

# =============================================================================
# Import Statements
# =============================================================================
"""
Pip install commands necessary for imported packages:
    * pip install numpy
    * pip install imageio
    * pip install imageio
    * pip install matplotlib
    * pip install tqdm
"""
import numpy as np
import argparse
import imageio
import logging
import matplotlib.pyplot as plt
from tqdm import trange

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

# =============================================================================
# Argument parser
# =============================================================================
def parse_args():
    '''
        * parse_args
        * Parses command line arguments for simulation parameters.
        * Equations:
            * Kinematic Relation (dimensionless):
                * dx'/dt' = v'
            * Dynamic Relation (dimensionless):
                * dv'/dt' = E'
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs:
                * None
            * Outputs:
                * args (argparse.Namespace): parsed command-line arguments

    The argument parser lets you override default simulation parameters
    via command-line flags, so you can easily adjust domain size, particle
    counts, time step, resolution, and other settings without editing the code.

    Default Initial Conditions:
        * Domain Length (L)                 = 1.0
        * Number of immobile ions           = 1000
        * Number of electrons               = 1000
        * Time step (dt)                    = 0.01
        * Total iterations (steps)          = 1000
        * Barnes-Hut opening angle (theta)  = 0.5
        * Grid resolution (N_grid)          = 200
        * Electron thermal velocity (v_th)  = 1.0
        * Random seed for reproducibility   = 0
        * Softening parameter (epsilon)     = 1e-3
        * Boundary decay constant (alpha)   = 200.0
        * Output animation filename         = "1D_Sheath_Treecode.mp4"
        * Snapshot indices                  = [0, steps//2, steps-1]
    '''
    parser = argparse.ArgumentParser(
        description="1D Plasma Sheath Simulation with Barnes–Hut Treecode"
    )
    parser.add_argument(
        "--version", action="version", version="plasma_sheath.py beta",
                        help="Show program version and exit.")
    parser.add_argument("--L", type=float, default=1.0,
                        help="Domain length (default: 1.0)")
    parser.add_argument("--N_ion", type=int, default=1000,
                        help="Number of immobile ions (default: 1000)")
    parser.add_argument("--N_electron", type=int, default=1000,
                        help="Starting number of electrons before boundary removal (default: 1000)")
    parser.add_argument("--dt", type=float, default=0.01,
                        help="Time step value (default: 0.01)")
    parser.add_argument("--steps", type=int, default=1000,
                        help="Total number of iterations (default: 1000)")
    parser.add_argument("--theta", type=float, default=0.5,
                        help="Barnes–Hut opening angle (default: 0.5)")
    parser.add_argument("--N_grid", type=int, default=200,
                        help="Grid resolution (default: 200)")
    parser.add_argument("--v_th", type=float, default=1.0,
                        help="Electron thermal velocity (default: 1.0)")
    parser.add_argument("--seed", type=int, default=0,
                        help="Random seed for reproducibility (default: 0)")
    parser.add_argument("--eps", type=float, default=1e-3,
                        help="Softening parameter to prevent singularity (default: 1e-3)")
    parser.add_argument("--alpha", type=float, default=200.0,
                        help="Decay constant for boundary correction (default: 200.0)")
    parser.add_argument("--output", type=str, default="1D_Sheath_Treecode.mp4",
                        help="Name of output MP4 file (default: 1D_Sheath_Treecode.mp4)")
    parser.add_argument("--snapshots", type=int, nargs = "+", default=None,
                        help=("List of time-step indices at which to save a static snapshot."
                              "If omitted, defaults to [0, steps//2, steps-1]."))
    args = parser.parse_args()
    # Validate
    if args.steps <= 0:
        parser.error("--steps must be positive.")
    if args.dt <= 0:
        parser.error("--dt must be positive.")
    if args.snapshots:
        invalid = [s for s in args.snapshots if s < 0 or s >= args.steps]
        if invalid:
            parser.error(f"Snapshot indices {invalid} out of range [0, {args.steps - 1}]")
    return args

# =============================================================================
# Treecode Implementation
# =============================================================================

def build_tree(xmin, xmax, pos, ch):
    '''
        * build_tree
        * Constructs a 1D Barnes–Hut tree for charge distribution to
        enable fast potential calculation.
        * This tree partitions space hierarchically and stores total
        charge and charge-weighted center.
        * Poisson's Equation (dimensionless):
            * ∂²φ'/∂x'² = -ρ'
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs:
                * xmin (float): left boundary of the current node region
                * xmax (float): right boundary of the current node region
                * pos (np.ndarray): particle positions within this node [1D array]
                * ch (np.ndarray): particle charges corresponding to pos [1D array]
            * Outputs:
                * Node (object): root of recursively built Barnes–Hut tree
    '''
    class Node:
        def __init__(self, xmin, xmax, positions, charges):
            self.x_min, self.x_max = xmin, xmax
            self.pos, self.ch = positions, charges
            self.q_total = np.sum(charges)
            self.x_cm = (
                np.sum(positions * charges) / self.q_total
                if self.q_total != 0 else 0.5 * (xmin + xmax)
            )
            self.left = self.right = None
            if len(positions) > 1 and not np.all(positions == positions[0]):
                midpoint = 0.5 * (xmin + xmax)
                left_mask  = positions <= midpoint
                right_mask = positions > midpoint
                if np.any(left_mask):
                    self.left  = Node(xmin, midpoint, positions[left_mask], charges[left_mask])
                if np.any(right_mask):
                    self.right = Node(midpoint, xmax,   positions[right_mask], charges[right_mask])
    return Node(xmin, xmax, pos, ch)

# =============================================================================
# Electrostatic Potential Energy Computation via Treecode Implementation
# =============================================================================
def compute_phi(x, node, theta, eps):
    '''
        * compute_phi
        * Recursively computes the electrostatic potential at position x using
        a Barnes–Hut tree with softening to prevent singularities.
        * Applies multipole approximation if the node is far away compared to its size.
        * Dirichlet BC Green's Function used separately.
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs:
                * x (float): position at which to compute the potential, evaluation point x'
                * node (Node): Barnes–Hut tree node
                * theta (float): opening angle threshold
                * eps (float): softening parameter
            * Outputs:
                * φ'(x') (float): computed potential at x
    '''
    if node is None:
        return 0.0

    dist = x - node.x_cm
    # Leaf or very close → direct sum
    if abs(dist) < eps or (node.left is None and node.right is None):
        return np.sum(node.ch / np.sqrt((x - node.pos)**2 + eps**2))

    dxn = node.x_max - node.x_min

    # Multipole approximation
    if dxn / abs(dist) < theta:
        return node.q_total / np.sqrt(dist**2 + eps**2)

    return (
        compute_phi(x, node.left,  theta, eps)
      + compute_phi(x, node.right, theta, eps)
    )

# =============================================================================
# Green's Function for Dirichlet BCs
# =============================================================================

def dirichlet_green_1d(x, xp, L):
    '''
        * dirichlet_green_1d
        * Computes the Green's function for a 1D domain [0, L] with Dirichlet
        boundary conditions: φ(0) = φ(L) = 0.
            * G(x',x'_p) for φ'(0)=φ'(1)=0
                * G = x'(1–x'_p) if x'<x'_p, else x'_p(1–x')
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs:
                * x (array_like): evaluation points x'
                * xp (array_like): source positions x'_p
                * L (float): domain length (=> 1 in primed)
            * Outputs:
                * G (np.ndarray): Green's function values at x
    '''
    xp_arr = np.array(xp)
    return np.where(
        x < xp_arr,
        x * (L - xp_arr) / L,
        xp_arr * (L - x)    / L
    )

# =============================================================================
# Boundary Correction Terms
# =============================================================================

def compute_AB(positions, charges, alpha=200.0, L=1.0):
    '''
        * compute_AB
        * Computes image charge corrections A and B to enforce grounded Dirichlet
            * A' = Σ q'_i G_left(x'_i), B' = Σ q'_i G_right(x'_i)
        boundary conditions in the domain via an exponential model.
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs:
                * positions (np.ndarray): particle positions [1D array], x_i'
                * charges (np.ndarray): corresponding particle charges [1D array], q_i'
                * alpha (float): decay constant for image-charge model
                * L (float): domain length
            * Outputs:
                * A' (float): left boundary correction coefficient
                * B' (float): right boundary correction coefficient
    '''
    xi     = 2 * positions / L - 1
    G_left = 0.5 * alpha * np.exp(-alpha * np.abs(-1 - xi))
    G_right= 0.5 * alpha * np.exp(-alpha * np.abs( 1 - xi))
    return np.sum(charges * G_left), np.sum(charges * G_right)

# =============================================================================
# Main Simulation Loop
# =============================================================================

def main():
    '''
        * main
        * Entry point: sets up simulation parameters, initializes particles,
        and runs the time-stepping loop to build trees, compute potentials,
        apply boundary corrections, and save animation.
            * for each step:
                1. build_tree → φ'_grid via compute_phi
                2. φ'_tot = φ'_grid - [(1–x')A' + x'B'] + sheath terms
                3. save snapshots if requested
                4. write frame to MP4
        * Date Created: 07/01/25
        * Date Last Modified: 07/21/25
        * Author: Beckett Henderson
        * API:
            * Inputs: args from parse_args
            * Outputs: animation + snapshot PNGs
    '''
    args = parse_args()
    np.random.seed(args.seed)

    # Snapshot indices not specified
    if args.snapshots is None:
        args.snapshots = [0, args.steps // 2, args.steps - 1]

    # Initial Particle Distributions
    ion_positions      = np.linspace(0, args.L, args.N_ion, endpoint=True) + args.L/(2*args.N_ion)
    ion_charges        = np.ones(args.N_ion)
    electron_positions = np.random.rand(args.N_electron) * args.L
    electron_charges   = -np.ones(args.N_electron)

    sheath_left_charge  = 0.0
    sheath_right_charge = 0.0

    # Prepare video writer
    writer = imageio.get_writer(args.output, fps=30)

    # Grid for potential & field
    grid_x = np.linspace(0, args.L, args.N_grid)

    # Precompute A, B once
    A, B = compute_AB(
        np.concatenate([ion_positions, electron_positions]),
        np.concatenate([ion_charges,   electron_charges]),
        alpha=args.alpha, L=args.L
    )

    # Time‑stepping loop
    for step in trange(args.steps, desc = "Simulating"):
        # Build tree for ions + electrons
        all_pos = np.concatenate([ion_positions, electron_positions])
        all_ch  = np.concatenate([ion_charges,   electron_charges])
        root    = build_tree(0.0, args.L, all_pos, all_ch)

        # Compute potential on grid
        phi_grid = np.array([
            compute_phi(x, root, args.theta, args.eps)
            for x in grid_x
        ])
        # Add image‐charge corrections & sheath
        phi_tot = (
            phi_grid
          - ((1 - grid_x/args.L)*A + (grid_x/args.L)*B)
          + sheath_left_charge  * dirichlet_green_1d(grid_x, 0.0,  args.L)
          + sheath_right_charge * dirichlet_green_1d(grid_x, args.L, args.L)
        )
        phi_tot[0] = phi_tot[-1] = 0.0

        # Snapshot Indices Specified -> PNG output
        if step in args.snapshots:
            plt.figure()
            plt.plot(grid_x, phi_tot)
            plt.title(f"φ at step {step}")
            plt.xlabel("x")
            plt.ylabel("φ")
            plt.savefig(f"snapshot_step{step}.png", dpi=150)
            plt.close()

    writer.close()
    print(f"Saved animation to {args.output}")

if __name__ == "__main__":
    main()