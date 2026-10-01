"""In this script, we read (for example) react.dump, whose atom
types are Si, O, H, Al. The interval of dump is (for example) 10 ps.
For every (for example) 100 snapshots, make neighbor lists,
Find all Al-O-Si or Al-O-Al, or Si-O-Si bonds,
and detect if the bonds are changed between the snapshots.
Write the result in (for example) bond_change.out :
    Snapshot   Time (ps)    changed Al-O-Si
    0          0            Al(12)-O(3)-Si(20) (O xyz: 1.0 2.0 3.0 Ang)     # here, 12 is 0-based index
    100        1000         Al(12)-O(3)  Si(20)  Si(22)-O(4)-Si(23) (O xyz: 4.0 5.0 6.0 Ang)
    200        2000         Si(22)-O(4)  Si(23) 

cutoff: Si-O: 2.172 Ang (Ovito GUI)
        Al-O: 2.112

A bond is broken only when two bonded atoms are apart by (for example) 1 Ang more than cutoff.
Also, a bond is formed only when two atoms were previously apart by (for example) 1 Ang (the same buffer as above line) more than cutoff.


Keep this doc because it's written by a human.
Written by codex, iterated with ovito gui.
"""

import argparse
from pathlib import Path

import numpy as np
from ase.neighborlist import neighbor_list

from thermo_SiO2.io import read_sil


def parse_arguments():
    """Return command-line settings for the bridge-change analysis."""
    parser = argparse.ArgumentParser(
        description="Find changes in Si/Al--O--Si/Al bridges in a LAMMPS dump."
    )
    parser.add_argument("--dump", default="react.dump", help="input LAMMPS dump")
    parser.add_argument(
        "--output", default="bond_change.out", help="output text file"
    )
    parser.add_argument(
        "--stride", type=int, default=100,
        # "--stride", type=int, default=1000,
        help="analyse every Nth dump snapshot (default: 100)",
    )
    parser.add_argument(
        "--si-o-cutoff", type=float, default=2.172,
        help="Si--O bond cutoff in Angstrom (default: 2.172)",
    )
    parser.add_argument(
        "--al-o-cutoff", type=float, default=2.112,
        help="Al--O bond cutoff in Angstrom (default: 2.112)",
    )
    parser.add_argument(
        # "--break-buffer", type=float, default=1.0,
        "--break-buffer", type=float, default=2.0,
        help="additional separation beyond the formation cutoff required to break a bond (default: 2.0 Angstrom)",
    )
    parser.add_argument(
        "--timestep-ps", type=float, default=0.001,
        help="MD timestep in ps per LAMMPS step (default: 0.001)",
    )
    parser.add_argument(
        "--symbols", default="Si O H Al",
        help="LAMMPS atom-type order passed to read_sil",
    )
    args = parser.parse_args()
    if args.stride < 1:
        parser.error("--stride must be at least 1")
    if args.si_o_cutoff <= 0 or args.al_o_cutoff <= 0:
        parser.error("--si-o-cutoff and --al-o-cutoff must be positive")
    if args.break_buffer < 0:
        parser.error("--break-buffer must be non-negative")
    if args.timestep_ps < 0:
        parser.error("--timestep-ps must be non-negative")
    return args


def bond_distances(atoms, si_o_cutoff, al_o_cutoff, break_buffer):
    """Return O--Si/Al distances that can satisfy either hysteresis threshold."""
    symbols = np.asarray(atoms.get_chemical_symbols())
    i_indices, j_indices, distances = neighbor_list(
        "ijd", atoms, max(si_o_cutoff, al_o_cutoff) + break_buffer
    )

    return {
        (int(i_atom), int(j_atom)): float(distance)
        for i_atom, j_atom, distance in zip(i_indices, j_indices, distances)
        if symbols[i_atom] == "O" and symbols[j_atom] in {"Si", "Al"}
    }


def bond_set(
    atoms, si_o_cutoff, al_o_cutoff, previous_bonds, previous_distances,
    break_buffer,
):
    """Return bonds using separate formation and breaking hysteresis rules."""
    distances = bond_distances(atoms, si_o_cutoff, al_o_cutoff, break_buffer)
    symbols = atoms.get_chemical_symbols()
    bonds = set()
    for bond, distance in distances.items():
        cutoff = {"Si": si_o_cutoff, "Al": al_o_cutoff}[symbols[bond[1]]]
        if bond in previous_bonds:
            if distance < cutoff + break_buffer:
                bonds.add(bond)
        elif (
            distance < cutoff
            and previous_distances.get(bond, float("inf")) > cutoff + break_buffer
        ):
            bonds.add(bond)
    return bonds, distances


def bridge_set(bonds):
    """Return bridge keys ``(O index, sorted cation-index pair)`` from bonds."""
    cations_by_oxygen = {}
    for oxygen, cation in bonds:
        cations_by_oxygen.setdefault(oxygen, set()).add(cation)

    bridges = set()
    for oxygen, cations in cations_by_oxygen.items():
        cations = sorted(cations)
        for first_position, first in enumerate(cations):
            for second in cations[first_position + 1:]:
                bridges.add((oxygen, first, second))
    return bridges


def format_bridge(bridge, symbols, positions):
    """Format a bridge key with its oxygen position for the output file."""
    oxygen, first, second = bridge
    first_symbol, second_symbol = symbols[first], symbols[second]
    left, right = sorted(((first_symbol, first), (second_symbol, second)))
    x, y, z = positions[oxygen]
    return (
        f"{left[0]}({left[1]})-O({oxygen})-{right[0]}({right[1]}) "
        f"(O xyz: {x:.6f} {y:.6f} {z:.6f} Ang)"
    )


def write_frame_events(handle, snapshot, atoms, formed, broken):
    """Write all changes detected for one sampled configuration."""
    timestep = atoms.info.get("timestep")
    if timestep is None:
        raise ValueError("trajectory frame has no LAMMPS timestep metadata")
    time_ps = timestep * write_frame_events.timestep_ps
    handle.write(
        f"Snapshot {snapshot}  Timestep {timestep}  Time (ps) {time_ps:.6f}\n"
    )
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    for event, bridges in (("formed", formed), ("broken", broken)):
        for bridge in sorted(bridges):
            handle.write(f"  {event:6s}  {format_bridge(bridge, symbols, positions)}\n")


def main():
    args = parse_arguments()
    dump_path = Path(args.dump)
    if not dump_path.is_file():
        raise FileNotFoundError(f"input dump does not exist: {dump_path}")

    cfgs = read_sil(
        str(dump_path), atom_symbols=args.symbols, index=f"::{args.stride}", silent=True
    )
    if not cfgs:
        raise ValueError(f"no configurations were read from {dump_path}")

    reference_symbols = cfgs[0].get_chemical_symbols()
    initial_distances = bond_distances(
        cfgs[0], args.si_o_cutoff, args.al_o_cutoff, args.break_buffer
    )
    previous_bonds = {
        bond for bond, distance in initial_distances.items()
        if distance < {"Si": args.si_o_cutoff, "Al": args.al_o_cutoff}[
            cfgs[0][bond[1]].symbol
        ]
    }
    previous_distances = initial_distances
    previous_bridges = bridge_set(previous_bonds)
    write_frame_events.timestep_ps = args.timestep_ps

    formed_count = broken_count = 0
    with Path(args.output).open("w", encoding="utf-8") as handle:
        handle.write("# Si/Al--O--Si/Al bridge changes\n")
        handle.write(
            f"# dump={dump_path} stride={args.stride} "
            f"Si-O_cutoff={args.si_o_cutoff:.6f} Ang "
            f"Al-O_cutoff={args.al_o_cutoff:.6f} Ang "
            f"break_buffer={args.break_buffer:.6f} Ang "
            "formation_requires_previous_separation_above_break_threshold=true "
            f"timestep_ps={args.timestep_ps:.9f}\n"
        )
        handle.write("# Atom indices are zero-based ASE indices.\n")
        for snapshot, atoms in zip(
            range(args.stride, len(cfgs) * args.stride, args.stride), cfgs[1:]
        ):
            if atoms.get_chemical_symbols() != reference_symbols:
                raise ValueError(
                    "atom symbols/order changed between frames; zero-based indices "
                    "cannot be compared safely"
                )
            current_bonds, current_distances = bond_set(
                atoms, args.si_o_cutoff, args.al_o_cutoff,
                previous_bonds, previous_distances, args.break_buffer,
            )
            current_bridges = bridge_set(current_bonds)
            formed = current_bridges - previous_bridges
            broken = previous_bridges - current_bridges
            if formed or broken:
                write_frame_events(handle, snapshot, atoms, formed, broken)
                formed_count += len(formed)
                broken_count += len(broken)
            previous_bridges = current_bridges
            previous_bonds = current_bonds
            previous_distances = current_distances

        handle.write(
            f"# analysed_frames={len(cfgs)} formed={formed_count} broken={broken_count}\n"
        )
    print(
        f"Analysed {len(cfgs)} frames; found {formed_count} formed and "
        f"{broken_count} broken bridges in {args.output}."
    )


if __name__ == "__main__":
    main()
