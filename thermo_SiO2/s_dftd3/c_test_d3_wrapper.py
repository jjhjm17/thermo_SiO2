"""This script tests the speed of d3 part for a test backbone structure."""

# ── Import ase.db BEFORE tensorflow/tensorpotential ──────────────────────────
from ase.io import read
from dftd3_wrapper import SDftd3Calculator
import time

if __name__ == '__main__':
    cfg = read('POSCAR')

    start = time.perf_counter()  # wall clock time
    n_repeats = 5
    # n_repeats = 1
    calc = SDftd3Calculator()  # Use the CLI wrapper instead of the direct DFTD3 calculator
    for i in range(n_repeats):
        cfg.positions[0, 0] += 0.001 * i  # Perturb the first atom slightly to avoid caching effects
        cfg.calc = calc

        print(f'{cfg.pbc =}')  # True for POSCAR
        calc.calculate(cfg, properties=['energy', 'forces'])
        energy = calc.results['energy']
        forces = calc.results['forces']
    end = time.perf_counter()
    print(f"Time taken for D3 calculation: {end - start:.3f} s")
    print(f"Average time: {(end - start) / n_repeats:.3f} s")
    print('Energy:', energy)
    print('Forces:', forces)
