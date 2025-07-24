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

- **Length scale:**                $L = \lambda_D$ (Debye length)  
- **Thermal velocity:**      $v_{th}$  
- **Time scale:**                  $T = L / v_{th}$  
- **Potential scale:**             $\phi_0 = \frac{m_e v_{th}^2}{q_0}$

With dimensionless variables:
- $x' = \frac{x}{L}$  
- $v' = \frac{v}{v_{th}}$  
- $t' = \frac{t}{T}$  
- $\phi' = \frac{\phi}{\phi_0}$  
- $E' = \left(\frac{q_0 L}{m_e v_{th}^2}\right) E$

---

## 🧪 Governing Equations

- **Kinematics:** $\frac{dx'}{dt'} = v'$
- **Dynamics:**  $\frac{dv'}{dt'} = E'$
- **Poisson Equation:**  $\frac{∂^2 \phi'}{∂x'^2} = -ρ'$
- **Leapfrog Updates:**
  - $v'_{n+1} = v'_n + E'(x') Δt'$
  - $x'_{n+1} = x'_n + v'_{n+1} Δt'$
- **Dirichlet Green's Function:**
  - $G(x', x'_p) = \begin{cases}
      x'(1 - x'_p), & \text{if } x' < x'_p \\
      x'_p(1 - x'), & \text{otherwise}
    \end{cases}$
- **Boundary Image Charge Terms:**
  - $A' = \sum q'_i G_{\text{left}}(x'_i)$  
  - $B' = \sum q'_i G_{\text{right}}(x'_i)$  

---

## ⚙️ Features

- Efficient **Barnes–Hut treecode** for 1D electrostatics
- Realistic **absorbing boundaries**
- **Command-line flags** to override all key parameters
- Generates an **MP4 animation** and optional **static snapshots**
- Reproducible results via **random seed control**

---

## 🚀 Getting Started

### 🔧 Install Dependencies

```bash
pip install numpy matplotlib imageio tqdm
```

### ▶️ Run the Simulation

```bash
python 1D_Sheath_Treecode.py
```
You can **customize** the simulation with optional arguments, for example:
```bash
python 1D_Sheath_Treecode.py --N_ion 500 --N_electron 800 --steps 1500 --output my_sim.mp4
```

### 📊 Example Output
- MP4:
- Snapshots:

### 📝 Default Parameters
| Parameter    | Value  | Description                    |
| ------------ | ------ | ------------------------------ |
| `L`          | 1.0    | Domain length                  |
| `N_ion`      | 1000   | Number of immobile ions        |
| `N_electron` | 1000   | Number of mobile electrons     |
| `dt`         | 0.01   | Time step                      |
| `steps`      | 1000   | Total time steps               |
| `theta`      | 0.5    | Barnes–Hut opening angle       |
| `N_grid`     | 200    | Number of grid points          |
| `v_th`       | 1.0    | Electron thermal velocity      |
| `eps`        | 1e-3   | Softening parameter            |
| `alpha`      | 200.0  | Boundary correction decay rate |
| `output`     | \*.mp4 | Animation file name            |

### 📁 Files Included
- 1D_Sheath_Treecode.py -- main simulation driver
- README.md -- project overview and usage instructions
- *.png -- optional output snapshots (generated during runtime)
- *.mp4 -- animation of sheath formation

## 📚 References
- **Christlieb, A.J., Krasny, R., & Verboncouer, J.P. (2004)**
  *Efficient particle Simulation of a virtual cathode using a grid-free treecode Poisson solver*
  *IEEE transactions on Plasma Science,* Vol. 32, No. 2
  🔗 ieeexplore.ieee.org/document/1308480
- **Christlieb, A.J., Krasny, R., Verboncouer, J.P., Emhoff, J.W., Boyd, I.D. (2006)**
  *Grid-free plasma Simulation techniques*
  *IEEE Transactions on Plasma Science,* Vol. 34, No. 2
  🔗 ieeexplore.ieee.org/document/1621283
- **Causley, M., Christlieb, A., Güçlü, Y., Wolf, E. (2013)**
  *Method of lines transpose: A fast implicit wave propagator*
- **Christlieb, Andrew J., Sands, William A., White, Stephen R. (2025)**
  *A Particle-in-Cell Method for Plasmas with a Generalized Momentum Formulation, Part I: Model   Formulation*
  *Journal of Scientific Computing,* Vol. 103, No. 15
  🔗 link.springer.com/article/10.1007/s10915-025-02824-1

## 💡Future Additions

## 🌐 License
This project is available for academic and educational use.
If you use this code in your research or work, please credit the author.
---
### ☑️ To use it:
1. Open a terminal in your repo folder.
2. Create the file:
   ```bash
   notepad README.md
   ```
3. Paste the content and save.
