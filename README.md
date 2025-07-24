# 1D Plasma Sheath Simulation with Barnes–Hut Treecode

Author: **Beckett Henderson**  
Date Created: **July 1, 2025**  
Last Modified: **July 21, 2025**

---

## 🌌 Overview

This project simulates a **1D plasma sheath** with **immobile ions** and **mobile electrons** using a **Barnes–Hut treecode** to efficiently compute electrostatic potentials. It applies **Dirichlet boundary conditions**, uses a **leapfrog integration scheme**, and visualizes sheath formation over time.

---

## 🧮 Non-Dimensionalization

The simulation is fully non-dimensionalized using the following reference scales:

- **Length:**                              $L = \lambda_D$ (Debye length)  
- **Thermal velocity:**                 $v_{th}$  
- **Time:**                             $T = L / v_{th}$  
- **Potential:**                      $\phi_0 = \frac{m_e v_{th}^2}{q_0}$

---

## 🧪 Governing Equations

- **Poisson’s Equation:**              $\frac{∂^2 \phi'}{∂x'^2} = -ρ'$
- **Leapfrog Integration:**
  - $v'_{n+1} = v'_n + E'(x') Δt'$
  - $x'_{n+1} = x'_n + v'_{n+1} Δt'$
- **Green’s Function (Dirichlet BC):**
  - $G(x', x'_p) = \begin{cases}
    x'(1 - x'_p), & \text{if } x' < x'_p \\
    x'_p(1 - x'), & \text{otherwise}
  \end{cases}$

---

## ⚙️ Features

- Efficient **Barnes–Hut treecode** for 1D electrostatics
- Realistic **absorbing boundaries**
- **Command-line flags** to override all key parameters
- Generates an **MP4 animation** and optional **static snapshots**
- Reproducible results via **random seed control**

---

## 🚀 Getting Started

### 🔧 Requirements

```bash
pip install numpy matplotlib imageio tqdm
