"""This script makes poreMS input, pore.yaml."""

import yaml
from math import pi
import numpy as np
from ase.io import read


def make_porems_yaml(dump_file, *, diameter):
    """This script makes poreMS input, pore.yml from a dump file.
    diameter in nm."""
    cfg = read(
        dump_file, format="lammps-dump-text", index=0, specorder="Si O H Al".split()
    )
    cell = cfg.get_cell() / 10  # nm
    # cell[i, j] is the jth Cartesian coordinate of the ith cell vector.
    # (vec_a + vec_b + vec_c) / 2
    centroid = (cell[0, :] + cell[1, :] + cell[2, :]) / 2
    length = float(np.linalg.norm(cell[2, :]))
    surface_area = float(2 * pi * diameter / 2 * length)
    rect_volume = float((2 * centroid)[0] * (2 * centroid)[1] * (2 * centroid)[2])
    # volume: volume of rectangular box
    dimensions = (2 * centroid).tolist()
    centroid = centroid.tolist()
    data = {
        "shape_00": {
            "diameter": diameter,
            "parameter": {
                "central": [0, 0, 1],
                "centroid": centroid,
                "diameter": diameter,
                "length": length,
            },
            "roughness": 0,
            "shape": "CYLINDER",
            "surface": surface_area,
            "volume": pi * (diameter / 2) ** 2 * length,
        },
        "system": {
            "centroid": centroid,
            "dimensions": dimensions,
            "reservoir": 0.000000001,
            "surface": {
                "ex": 0,
                "in": surface_area,
            },
            "volume": rect_volume,
        },
    }
    with open("pore.yml", "w") as outfile:
        yaml.dump(data, outfile, default_flow_style=False)


def write_lammps_dump_text(file, configs, time_step, box_bound_text):
    with open(file, "w") as fout:
        for index, config in enumerate(configs):
            fout.write("ITEM: TIMESTEP\n")
            fout.write(f"{index * time_step}\n")
            fout.write("ITEM: NUMBER OF ATOMS\n")
            fout.write(f"{len(config)}\n")
            fout.write(f"{box_bound_text}")
            fout.write("ITEM: ATOMS id type xu yu zu\n")
            if config.get_chemical_formula(mode="hill", empirical=True) != "O":
                print("Error: now only O atoms work. Please code more")
                sys.exit()
            for i_atom, atom in enumerate(config):
                pos = atom.position
                fout.write(f"{i_atom} 2 {pos[0]} {pos[1]} {pos[2]}\n")


def cell_hex_to_rect(dump_file, output_dump, time_step, box_bound_text):
    configs = read(
        dump_file,
        index=":",
        format="lammps-dump-text",
        specorder=["Si", "O", "H", "Al"],
    )
    # cell = configs[0].get_cell() / 10  # nm
    # # cell[i, j] is the jth Cartesian coordinate of the ith cell vector.
    # # (vec_a + vec_b + vec_c) / 2
    # # We assume hexagonal cell with and assume that
    # # the box_point can be found by this way.
    # centroid = (cell[0, :] + cell[1, :] + cell[2, :])/2
    # box_point = 2 * centroid
    #     box_bound_text = f"""ITEM: BOX BOUNDS xy xz yz pp pp pp
    # -1.1122981778246457e-01 7.6572057490665628e+01 0
    # -9.6327847857328050e-02 4.4176788722061595e+01 0.0000000000000000e+00
    # -5.5614908891232287e-02 2.5505480860591746e+01 0.0000000000000000e+00
    # """
    # manually change xy to 0

    # cell = configs[0].get_cell()
    # # print(f'{rectangle_cell = }')

    # for config in configs:
    #     config.set_cell(rectangle_cell)

    # print(f'{configs[0].get_cell() = }')
    write_lammps_dump_text(output_dump, configs, time_step, box_bound_text)
