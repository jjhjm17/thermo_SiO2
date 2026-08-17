"""
This script obtains IR spectrum of each OH groups.
We use formal charge of H: +1 to get the dipole moment from the OH vector,
which depends on the difference of OH positions.
Input
    in.yaml:
        OH_vector_H_indices : [700, 682]  # starts from 0
        OH_output_plot = "IR_spectrum_OH.pdf"
        OH_output_data = "IR_spectrum_OH.dat"
    OH_vectors.dat
Output
    OH_output_plot: a plot of all partial IR's of OH vectors,
        with legens indicating the index of H atoms,
        similar to output_plot of ir.py.
    OH_output_data: data file, similar to the output_data of ir.py,
        but for many intensities for each H in OH_vector_H_indices.

Details 
    After running get_OH_vectors.py, we have OH_vectors.dat
        ) head OH_vectors.dat
        # step H_700_x H_700_y H_700_z H_682_x H_682_y H_682_z  (OH vector components, Ang)
        0 0.57550 0.73990 0.22780 -0.93727 0.17790 -0.08117
        1 0.57575 0.75050 0.23230 -0.93041 0.18010 -0.08376
        ...


