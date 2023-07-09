# GAIT: A New Fingerprint

Welcome to our project repository! This project leverages the actibelt step detection algorithm and enhances it to calculate a 'norm step'. A norm step represents the result of averaging multiple steps, creating a unique identifier similar to a fingerprint. All codes in this project are written in the Julia programming language.

## Repository Structure

- `Plots/`: Contains all plots derived from applying our algorithm.
- `Our data/`: Includes our own data recordings used for the project.
- `GaitAnalysis.jl`: This is the final script that generates our final plots for our presentation. It is independent and does not require any other script to function. It operates on our data and the data from the actibelt database to create norm steps after filtering. Then, it calculates and visualizes their L2 distances using a heat map.
- `GaitAnalysis_previous_version.jl`: A previous version of the final script with variations in the plotting form and the filtering steps.
- `Test_aydin.jl`: A proto-code for the implementation of the norm step algorithm by Aydin Uzun. This script is included for documentation purposes only and is not necessary to run the code.

Please note, some data from the actibelt database was used during the project, but it is not included in this repository due to copyright reasons.

## Contact

If you have any questions or would like to discuss this project further, please feel free to reach out.

**Uzun, Aydin**  
Email: aydin.uzun@tum.de
