# Targeting SARS-CoV-2 Main Protease with D-peptides
Welcome to D-peptide binder design project 😃 This repository is the source code for Targeting SARS-CoV-2 Main Protease with D-peptides.
![workflow](https://github.com/laiyii/D-peptide-binder-design/blob/main/Dpep_fig1.png)
If you have any questions, feel free to discuss in [Issues](https://github.com/laiyii/D-peptide-binder-design/issues) or contact Laiyi (fly_ccme@pku.edu.cn).

## Installation
> **Note:** Ubuntu 20.04 is recommended for the procedure.
```shell
git clone https://github.com/laiyii/D-peptide-binder-design.git .
vim ~/.bashrc
export DPEP="/path/to/D-peptide-binder-design/source_code"
source ~/.bashrc
```
### Other applications in the workflow
- [Naccess](http://www.bioinf.manchester.ac.uk/naccess/)
- [Rosetta3.11](https://downloads.rosettacommons.org/software/academic/)
- [gmx_MMPBSA](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00645)
- [PASTA2.0](https://doi.org/10.1093/nar/gku399)

## Tutorial
### Prepare curled L-helical scaffolds
You can generate scaffolds with customized needs.
```shell
gcc $DPEP/curled_lib/script/PhiPsi2Helix.c -o $DPEP/curled_lib/script/PhiPsi2Helix -lm
```
Running PhiPsi2Helix to generate scaffolds with given parameters:
```shell
chmod +x $DPEP/curled_lib/script/curl_helix_gen.sh
$DPEP/curled_lib/script/curl_helix_gen.sh -o H_-62_-39_-3_-1_-60.pdb -outdir helix_lib -len 21 -phi0 -62.0 -delphi -3.0 -psi0 -39.0 -delpsi -1.0 -phase -60.0 --
```
where `-o` and `-outdir` defines the output name and output directory, `-len` is the length of the polyALA sequence. Rational range of other parameters are shown in [Table S1].<br>
We also provide helix scaffold library with various lengths (21 aa, 24 aa, 28 aa, 31 aa, 35 aa, 38 aa, and 42 aa) already generated in this work. Click [here](https://1drv.ms/u/c/1838b20033e25fae/EcgmP7MWDtxGiOSvWAjSSzwBrgVcsVyyKKK8k4YAJU5nkg?e=Xm1Qxd) to download.

### Docking of the helical scaffolds to the target
#### Flip the target into D-type
Before docking, please flip your target to D-type, with residue names unchanged. Note that input file type should be a pdb file **with hydrogens removed**.
```shell
chmod +x $DPEP/docking/mirror_target/mirror_target.sh
$DPEP/docking/mirror_target/mirror_target.sh -i your_input_file.pdb -o your_output_file.pdb
```
The default output of `-o` is your_input_file_mirror.pdb
#### Surface residues remark
To generate grid scores, we need to define surface atoms (with atom-wise SASA larger than 1 Å²).
```shell
./naccess your_input_file_mirror.pdb
```
> **Note:** For usage of naccess, see [Naccess homepage](http://www.bioinf.manchester.ac.uk/naccess/).

With the asa file, we will edit the target structure file at column 70. Surface atoms are marked as 1, while internal atoms are 0: 
```shell
python3 $DPEP/docking/mirror_target/surf_protein.py -i_asa mono_mirror_noh.asa -i your_input_file_mirror.pdb -o your_input_file_mirror_surf.pdb
```
`-i_asa` is the output file of Naccess.

#### Helix scaffolds docking
Compile the file first.
```shell
gcc $DPEP/docking/HelixScaffoldDocking/HelixScaffoldDocking_batch.c -o $DPEP/docking/HelixScaffoldDocking/HelixScaffoldDocking_batch -lm
```
Running the program with:
```shell
chmod +x ./$DPEP/docking/HelixScaffoldDocking/HSD_batch.sh
./$DPEP/docking/HelixScaffoldDocking/HSD_batch.sh -t target_processed.pdb -b batch_info -a central_atom_id
```
Docking tasks are performed in batch. The input options include:<br>
- `-t` Processed target structure.
- `-b` A text file containing input scaffolds location and output file names.
- `-a` Atom ID, the ID of atom as the center of docking box.
<br>
Here's an example for batch_info:

```text
$DPEP/docking/HelixScaffoldDocking/batch_info_example
```

The input scaffold file and output file is separated by spaces.

### Loop modeling with CCD
The following adjustments need to be made to the output structure.
1. Replace the target structure with the initial structure (containing H atoms) and name it to chain B.
2. Add H atoms to polyALA scaffold (ligand) and name it to chain A.
3. Rearrange the complex structure, put ligand to the first.
Then running loop modeling:
```shell
loopmodel.mpi.linuxgccrelease @ccd.flags
```
You can edit settings based on your task. Input files for CCD loop modeling is here:<br>
```text
$DPEP/docking/loop_modeling
```
### Sequence design
Before sequence design, residue name of target (chain B) should be changed.
```shell
python3 $DPEP/sequence_design/change_chainB_to_D_type.py -i input_file.pdb
```
Sequence design is excecuted by RosettaScripts. A template xml file is in:<br>
```text
rosetta_scripts.mpi.linuxgccrelease @$DPEP/sequence_design/Dpep_design.flags
```

### Sequence selection
Criteria for *in silico* sequence selection is described in the [paper](). Binding energy is calculated by gmx_MMPBSA and free energy of aggregation is calculated by PASTA2.0.<br>
> **Note:** Complex structure is flipped to D-ligand and L-receptor after sequence design.

## About database
There are two csv files in the branch. <br>
`geometry_score_database.csv` is the database used to collect geometry characters of helical ligands, PDB IDs are listed in the first column.<br>
`interface_propensity_score_database.csv` is the database used to construct interface propensity scores. The first column is the PDB id, each structure is split into two parts to analyze ratios of different atoms at interface or surface, represented as chain IDs in the second column.

## License

This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this project in your research, please cite it using the following BibTeX entry:

```bibtex
@misc{your_project_name,
  author = {Your Name},
  title = {Your Project Title},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/yourusername/your-repository}},
}
```

## Acknowledgements
This design procedure involves multiple softwares related to protein design. We acknowledge and thank the developers of Naccess, Rosetta, gromacs, gmx_MMPBSA and PASTA2.0 for their incredible and hard work.

