# Targeting SARS-CoV-2 Main Protease with D-peptides
Welcome to D-peptide binder design project 😃 This repository is the source code for Targeting SARS-CoV-2 Main Protease with D-peptides.
![workflow](https://github.com/laiyii/D-peptide-binder-design/blob/main/Dpep_fig1.png)
If you have any questions, feel free to discuss in [Issues](https://github.com/laiyii/D-peptide-binder-design/issues) or contact Laiyi (fly_ccme@pku.edu.cn).

## Installation
```shell

```
### Other applications in the workflow
- [Naccess](http://www.bioinf.manchester.ac.uk/naccess/)
- [Rosetta3.11](https://downloads.rosettacommons.org/software/academic/)

## Tutorial
### Prepare curled L-helical scaffolds
You can generate scaffolds with customized needs.
```shell
gcc $HSD/curled_lib/script/PhiPsi2Helix.c -o $HSD/curled_lib/script/PhiPsi2Helix -lm
```
Running PhiPsi2Helix to generate scaffolds with given parameters:
```shell
chmod +x $HSD/curled_lib/script/curl_helix_gen.sh
$HSD/curled_lib/script/curl_helix_gen.sh -o H_-62_-39_-3_-1_-60.pdb -outdir helix_lib -len 21 -phi0 -62.0 -delphi -3.0 -psi0 -39.0 -delpsi -1.0 -phase -60.0 --
```
where `-o` and `-outdir` defines the output name and output directory, `-len` is the length of the polyALA sequence. Rational range of other parameters are shown in [Table S1].<br>
We also provide helix scaffold library with various lengths (21 aa, 24 aa, 28 aa, 31 aa, 35 aa, 38 aa, and 42 aa) already generated in this work. Click [here](https://1drv.ms/u/c/1838b20033e25fae/EcgmP7MWDtxGiOSvWAjSSzwBrgVcsVyyKKK8k4YAJU5nkg?e=Xm1Qxd) to download.

### Docking of the helical scaffolds to the target
#### Flip the target into D-type
Before docking, please flip your target to D-type, with residue names unchanged. Note that input file type should be a pdb file **with hydrogens removed**.
```shell
chmod +x $HSD/docking/mirror_target/mirror_target.sh
$HSD/docking/mirror_target/mirror_target.sh -i your_input_file.pdb -o your_output_file.pdb
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
python3 $HSD/docking/mirror_target/surf_protein.py -i_asa mono_mirror_noh.asa -i your_input_file_mirror.pdb -o your_input_file_mirror_surf.pdb
```
`-i_asa` is the output file of Naccess.

#### Helix scaffolds docking
Compile the file first.
```shell
gcc $HSD/docking/HelixScaffoldDocking/HelixScaffoldDocking_batch.c -o $HSD/docking/HelixScaffoldDocking/HelixScaffoldDocking_batch -lm
```
Running the program with:
```shell
chmod +x ./$HSD/docking/HelixScaffoldDocking/HSD_batch.sh
./$HSD/docking/HelixScaffoldDocking/HSD_batch.sh -t target_processed.pdb -b batch_info -a central_atom_id
```
Docking tasks are performed in batch. The input options include:<br>
- `-t` Processed target structure.
- `-b` A text file containing input scaffolds location and output file names.
- `-a` Atom ID, the ID of atom as the center of docking box.
Here's an example for batch_info:
```text
$HSD/docking/HelixScaffoldDocking/batch_info_example
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
Required files for CCD loop modeling is in:<br>
```text
$HSD/docking/loop_modeling
```
### Sequence design
Before sequence design, residue name of target (chain B) should be changed.
```shell
python3 $HSD/sequence_design/change_chainB_to_D_type.py -i input_file.pdb
```
Sequence design is excecuted by RosettaScripts. A template xml file is in:<br>
```text
rosetta_scripts.mpi.linuxgccrelease @$HSD/sequence_design/Dpep_design.flags
```

## License

## Acknowledgements






