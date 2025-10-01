If you use this code, please cite:

Fu F., Zhuo R., Chen X. Social Imitation Dynamics of Vaccination Driven by Vaccine Effectiveness and Beliefs. 

# Lattice Vaccination–Epidemic Simulations (C++)

Simulations of coupled behavior–disease dynamics on 2D lattices. Each season consists of:

1. **Stage 1 (Behavior):** Individuals choose whether to vaccinate. Social imitation uses a Fermi update with intensity \(K\) (bounded rationality).
2. **Stage 2 (Epidemic):** SIR dynamics on the lattice. Vaccination reduces susceptibility by factor \((1-\varepsilon)\). Recovery rate \(\gamma\), transmission rate \(\beta\), \(I_0\) initial seeds.

The code generates CSV outputs for vaccination coverage, epidemic size, and hysteresis analyses as vaccine effectiveness \(\varepsilon\) or perceived costs vary.

---

## Features
- Square lattice \(L \times L\) with periodic boundary conditions.
- Seasonal two-stage dynamics (behavior → epidemic) repeated until convergence.
- Imitation via Fermi rule with intensity \(K\).
- Vaccine effectiveness \(\varepsilon\) (primary failure).
- Fixed fraction of vaccine-averse agents.
- CSV outputs ready for plotting (Python/R/Matlab).

---

## Quick Start

### Requirements
- C++17 compiler (g++ 9+, clang 10+, MSVC 2019+)
- CMake ≥ 3.15
- (Optional) OpenMP for multithreading

### Build
```bash
git clone https://github.com/<org>/<repo>.git
cd <repo>
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
