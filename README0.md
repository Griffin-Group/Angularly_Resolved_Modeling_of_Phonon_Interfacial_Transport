# Example Project

This is a repository for my project on the  properties of actinide pyrochlores under strain. The preprint is available on arXiv (arxiv:7212.12345).

## Computational Info
All work in this repository was performed by me, except for the slab calculations (Spike Spiegel, Griffin Group), and the calculation of the MCAE (Dr. Jet Black, Bebop Group, Mars Autonomous University). All calculations were performed on the Swordfish II cluster at the University of Ganymede. We used VASP version 11.7.5 with the SecretNuclearWeapons patch (version 0.5) available for free on GitHub. We found that using more than 25 thousand GPUs with these calculations leads to a variety of MPI-related crashes, but the parallelism settings in the job script and `INCAR` should work fine.

## Structures
For Bk2Sn2O7, I carefully relaxed the experimental bulk structure (from ICSD-247234) and used it in all subsequent bulk calculations. For Cf2Ti2O7, I used the structure from the Materials Project (mp-43242432, retrieved December 2071) without further relaxation. The slabs were constructed by forming a 9-layer slab and fully relaxing it, then sandwiching in 2 extra layers near the center and selectively relaxing only the surrounding atoms (to save computational costs).

## Provenance
All non-self-consistent calculations started from the charge density of the self-consistent calculation (`1.scf`) in the same directory. The Bk2Sn2O7 calculations were easy to converge by starting from a colinear calculation (included in the `1a.scf_col` folder) but the Cf2Ti2O7 calculations were quite difficult and expensive due to the geometrical frustration. I backed up the `WAVECAR` and `CHGCAR` files to HPSS in the `Dy2Ti2O7/bulk/1.relax` directory. 

## Scripts
The Julia script to calculate Z54 topological invariants (`scripts/Z54_calc.jl`) was provided by Dr. Asimov Solensan (Jupiter Institue of Theoretical Physics) and we were asked not to share it outside the group without their permission.

## Other Files
The custom `POTCAR`s for Bk and Am are available in group directories. They were purchased from a local arms dealer on Callisto, and violate export controls, but the other `POTCAR`s are standard. 

## Figures
The band structure plots in Fig. 2 were generated using `sumo` version 20.7 with a small color patch that will be merged upstream soon (see the pull request here: [here](?))

---

# Data Policy

This is a sample repository for Griffin Group papers/projects. Group members are expected to create such repositories at the end of a project or before leaving the group for the following reasons:

1. **Knowledge Transfer**. Knowing converged parameters, expected results etc. saves other group members precious time, gives us something to compare our results to, and somewhere to start from.
2. **Reproducibility**. Being able to reproduce former group members' calculations with minimal effort after their departure is critical, especially if taking over someone's paper that's under review.

We want our science to be easily reproducible, and we don't want to reinvent the wheel. To this end, the expectation is that you create a repository for your projects with the information and files necessary to meet those two goals. Constructing such a repository should be easy, and we leave the structure up to you. You will mainly need to:

1. Write a brief `README.md`
2. Collect and organize the required calculation files. Many group members already use directory trees similar to what we suggest below.
3. Collect your analysis scripts and notebooks. You don't need to clean them up, just make sure they're readable and have minimal comments.
4. If you're leaving the group before your manuscript is accepted for publication: include your figure files and plotting scripts/notebooks.

# Recommended Directory Structure
We recommend you structure your repository in a "self-documenting" way to minimize the amount of work you have to put into writing the `README`. This example repository is structured as follows:

```
.
├── Cf2Ti2O7
│   ├── equil
│   │   ├── 0.relax
│   │   ├── 1.scf
│   │   ├── 2.bands
│   │   └── 3.dos
│   ├── strain_+0.5
│   │   ├── 0.relax
│   │   ├── 1.scf
│   │   ├── 2.bands
│   │   └── 3.dos
│   └── strain_-0.5
│       ├── 0.relax
│       ├── 1.scf
│       ├── 2.bands
│       └── 3.dos
├── Bk2Sn2O7
│   ├── 001_surface
│   │   ├── 1.scf
│   │   └── 2.bands
│   └── bulk
│       ├── 0.relax
│       ├── 1a.scf_col
│       ├── 1b.scf
│       ├── 2.bands
│       └── 3.dos
├── figures
│   ├── fig1.ai
│   ├── fig1_a.png
│   ├── fig1_b.png
│   ├── fig2.ai
│   ├── fig2_a.eps
│   ├── fig2_b.eps
│   ├── fig2_c.eps
│   ├── fig3.ai
│   └── fig3.pdf
├── scripts
│    ├── MCAE_analysis.ipynb
│    ├── Z54_calc.jl
│    ├── plot_bulk.py
│    └── plot_surface.py
└── README.md
```

At the top level, you have a directory for each material studied in the project. If you studied multiple structures of a given material (For example: bulk and several surfaces, conventional and primitive bulk, the same bulk cell under different strains, etc.), create directories for each structure. Within those directories, create numbered directories for each full calculation. The numbers should help naturally establish provenance. For example, the `2.bands` directory can be safely assumed to contain a non-self-consistent band structure calculation that used the charge density of the self-consistent calculation in `1.scf`, and so on. It's up to you how you do this. See below for the required files for each calculation.

You should also include a directory for analysis/plotting scripts/notebooks, a directory for figure files, and a top-level `README.md` file. See below for guidance on which files to include.

# `README.md`
You should include a top-level `README.md` in your repository (such as the first section of this file) that provides any information that cannot be directly gleaned from the files in the repository. This includes:

- [ ] A sentence or two about computational setup (version, any hacks/patches you used, which cluster, etc.)
- [ ] Where do your structures come from? 
    - [ ] If you did not relax your structures, where did you get them from?
    - [ ] If they were relaxed, it is still sometimes helpful to know where you started (i.e., ICSD ID or Materials Project ID).
    - [ ] If you include an MP ID, please include an approximate retrieval date since the structures get updated.
    - [ ] If you built slabs/nanowires/etc. from a bulk structure, or if you did any fancy multi-step relaxations (see README example), add a sentence or two about it.
- [ ] Explain the provenance of your calculations. 
    - [ ] Where did the charge density of any non-self-consistent calculations come from? If you used the recommended directory structure, this should be self-evident, but it's still helpful to mention it.
    - [ ] If you started a self-consistent calculation from a wave function/charge density to help convergence (e.g., starting a non-colinear calculation from a colinear `WAVECAR`), please mention it.
- [ ] If there are any relevant files on HPSS, please mention their location.
- [ ] If certain aspects of the project were performed by someone else (e.g., a collaborator or another group member), please mention it.
- [ ] Include any other information worth knowing for someone trying to reproduce the calculations or taking over the project. For example, did you find any weird behavior in certain VASP versions? Any odd bugs with certain settings?
- [ ] If there is no manuscript yet (e.g., you are graduating before a project is finished), please include any details that would normally go into a manuscript.

# Calculation Files
These are the files for calculations (e.g., DFT, GW, etc.) that we ask everyone to include in their repositories.

## All Codes
No matter what code you use, please include the following for each calculation:
- [ ] **Required:** `stdout`/`stderr` file. 
- [ ] **Required:** job script. If the calculation was run interactively, please include your `srun`/`mpirun` or equivalent command with necessary environment variables in a shell script.

Some codes (e.g., Quantum ESPRESSO) write their full output to `stdout`, and others (e.g., `VASP`) provide useful information and warnings. While the main parallelization options in a job script can be reconstructed from some output files, providing the file will also tell us what environment variables you used (some can affect performance or whether a calculation runs at all) and the exact executable you ran (i.e., was it a specific module version or was it a binary compiled with a certain compiler?)

If you are using a custom binary, please include the path in the job script instead of your `.bashrc`.

If your calculation hit a "hidden" bug, knowing exactly how you ran it is critical for reproducibility. (For example, VASP <6.4 has a hard-to-track bug in magnetic symmetry when compiled with PGI, but not with Intel.)

## VASP
- [ ] **Required**: `INCAR`
- [ ] **Required**: `KPOINTS`
- [ ] **Required**: `POSCAR`
- [ ] **Required**: `POTCAR.spec` (don't include the `POTCAR`)
- [ ] **Required**: `OUTCAR`
- [ ] **Possibly Required**: `WAVECAR` (on HPSS)
- [ ] **Possibly Required**: `CHGCAR` (on HPSS)

POTCAR files are proprietary and they can't be shared. Instead, we can generate a `POTCAR.spec` file:

```bash
sed -n 's/^ *TITEL *=//p' POTCAR > POTCAR.spec
```
which will generate something that looks like this:
```
  PAW_PBE Dy_3 06Sep2000                 
  PAW_PBE Ti_pv 07Sep2000                
  PAW_PBE O 08Apr2002                    
```
which tells us the exact `POTCAR` you used as well as the *version* of the `POTCAR`. (New VASP releases will sometimes improve pre-existing `POTCARs` to fix bugs or ghost-states.) The idea is similar to MP `POTCAR.spec` files, but those don't include the version/generation date.

`OUTCAR` files include much useful information but they're still lightweight and easy to share on GitHub. If your calculations were extremely hard to converge, you should back up the `WAVECAR` and `CHGCAR` to the group directory on HPSS. (Instructions coming soon.)

## Quantum ESPRESSO
- [ ] **Required**: input file (often named `pw.in` or `scf.pwi`, etc.)
- [ ] **Required**: output file (i.e., `stdout`/`stderr` file, already mentioned above.)
- [ ] **Required**: XML file (if under 50 MB)
- [ ] **Required**: UPF file for every pseudopotential used.

Quantum ESPRESSO's pseudopotentials are open source and can come from many libraries. It is recommended, but not required, that you keep the default name assigned by the library (e.g., `Pb.rel-pbe-dn-rrkjus_psl.1.0.0.UPF`).

## Other Codes
If you're using another DFT or post-processing code that is not listed here, we trust you to use your judgment and the guidelines outlined in this note to include the necessary files. If you had to rerun the calculation yourself in a year or two, what files would you need to do so? That's usually one or more input files and possibly some kind of data set (e.g., a pseudopotential or PAW set) or the output of another calculation that might be hard to reproduce (e.g., a `WAVECAR` from a very expensive calculation). If other lightweight files are helpful but not strictly necessary, you can also include them.

# Scripts

- [ ] **Required**: relevant scripts/notebooks.

If your project involved any kind of analysis in a script (beyond basic plotting or extracting numbers from output files), whether home-brew or using a well-established package, please include the script or notebook. There is no need to fully polish it, these repositories are intended to be private, but please make sure it's readable and has some basic comments.

If you used a script to only extract data from an output file, but you wrote the parser yourself because no open-source options were available, please include the script.

# Figures

You only need to worry about this if you are leaving the group **before the paper is accepted for publication**. In this case, another group member might have to deal with the peer-review process after your departure, and we want to minimize the work they have to put in.

- [ ] **Required**: Plotting scripts and notebooks (under `./scripts`)
- [ ] **Required**: Illustrator/InkScape (or similar) figure files (in `.ai`/`.svg` or similar)
- [ ] **Required**: Subfigure files used in the Illustrator files (in `.pdf`/`.eps` if possible, otherwise `.png` or similar)

The plotting scripts/notebooks should include a plotting example for each kind of figure in the paper (e.g., if your paper has 6 similar band structures, you only need to include a plotting example for one of them). If you used `sumo` or some other command-line plotting tool, please include the command you used to generate the figure in a shell script to ensure consistency.

For plots, it is strongly recommended that you do most standard modifications (e.g., adjusting line widths and colors, tick lengths and widths, etc.) in the plotting script itself, and only use Illustrator/InkScape for things like annotations, combining subfigures, adding labels, as well as more complex modifications. This makes it much easier to generate consistent figures if reviewers ask for modifications or new calculations to be added.

# Preparing your files
Let's say you have a directory `MessyCalcs` that looks somewhat like this repository (i.e., it has a bunch of material directories, each with subdirectories for each structure, each with multiple calculation directories). You can create a `POTCAR.spec` in every single folder with the following command:

```bash
find MessyCalcs -name 'POTCAR' -execdir bash -c "sed -n 's/^ *TITEL *=//p' '{}' > POTCAR.spec"  \;
```

Then to get rid of all the files you don't need, you can use `rsync` with include and exclude filters to create a new directory `CleanCalcs`

```bash
rsync -avr MessyCalcs/ CleanCalcs/ --include="*/" --include={INCAR,KPOINTS,POSCAR,POTCAR.spec,job.out,job.sh,OUTCAR} --exclude='*'
```

Note that `job.sh` and `job.out` is what I call all my shell scripts and `stdout`/`stderr` files, respectively. Modify these as needed. If your directory structure wasn't well organized to start with, you should now go clean it all up before creating the git repo and adding `README.md` and other files.
