# Capacitor Matching Algorithm 

This algorithm utilizes Constraint Programming via Google OR-Tools to solve a specific electrical engineering challenge: finding the optimal combination of four capacitors from a discrete inventory to reach a precise target capacitance.

It does not just look for numerical accuracy; it optimizes for network balance, minimizing the "mismatch" between components to ensure circuit stability and performance.

### The Network

The algorithm assumes a specific bridge-like configuration where capacitors are arranged in two parallel-series branches. This layout is common in RF matching networks or differential sensing applications.

The equivalent capacitance $C_{eq}$ is calculated based on the following relationship:

$$C_{eq} = \frac{(C_a + C_c) \cdot (C_b + C_d)}{(C_a + C_c) + (C_b + C_d)}$$

Where:

Slot A & C are in parallel.

Slot B & D are in parallel.

The two resulting pairs are then placed in series.

### Key Features

Discrete Solver Integration: Since the cp_model operates on integers, the algorithm automatically scales pF values and uses a "guard scale factor" to maintain high precision during division operations.

Multi-Objective Optimization: The solver minimizes a complex cost function that balances:

Hard Range Compliance: Heavily penalizes any solution falling outside the user-defined ERROR range.

Absolute Accuracy: Minimizes the delta between the actual equivalent capacitance and the TARGET.

Set Spread (Mismatch): Encourages $C_a$ to match $C_c$, and $C_b$ to match $C_d$ for better thermal and electrical symmetry.

Directional Balancing: Applies penalties to prevent asymmetrical configurations that could negatively impact circuit performance.

Detailed Reporting: Generates a CLI report showing component assignments for Slots A, B, C, and D, as well as percentage mismatches and final error.


### How it Works

1. Scaling Logic

To handle decimals within an integer-based solver, the script applies a scale factor. It also uses a limit to prevent integer overflow during multiplication steps in the capacitance formula.

2. Penalty System

The algorithm uses weighted priorities to guide the solver

3. Solving

The CapMatch.solve() method invokes the CP-SAT solver, which explores all possible combinations of the provided inventory (allowing for duplicates) until it finds the configuration with the lowest total penalty.

### Requirements

Python 3.x

ortools (pip install ortools)

numpy (pip install numpy)

### Usage

Simply run the script:

python match.py
