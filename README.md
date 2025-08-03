# Docking tools and scaffold libs for *Targeting SARS-CoV-2 Receptor Binding Domain and Main Protease with D-peptides*
This repository is the source code and tutorial for *Targeting SARS-CoV-2 Receptor Binding Domain and Main Protease with D-peptides*, including curved helical scaffold library generation and D-peptide docking program. If you have any questions, feel free to discuss in [Issues](https://github.com/laiyii/D-peptide-binder-design/issues).<br>
![workflow](https://github.com/laiyii/D-peptide-binder-design/blob/main/figs/Dpep_fig1.jpeg)


## Installation
> **Note:** Ubuntu 20.04 is recommended.
```shell
git clone https://github.com/laiyii/D-peptide-binder-design.git .
vim ~/.bashrc
export DPEP="/path/to/D-peptide-binder-design/source_code"
source ~/.bashrc
```

## Tutorial
### Curved helical scaffold library generation
You can generate scaffolds with customized needs.
```shell
gcc $DPEP/curled_lib/script/PhiPsi2Helix.c -o $DPEP/curled_lib/script/PhiPsi2Helix -lm
```
Running PhiPsi2Helix to generate scaffolds with given parameters:
```shell
chmod +x $DPEP/curled_lib/script/curl_helix_gen.sh
$DPEP/curled_lib/script/curl_helix_gen.sh -outdir <output_directory> -len <length> -paramlist <csv_file>
```
where `-outdir` defines the output directory, `-len` is the length of the polyALA sequence. Range of other parameters are defined in <csv_file> (see $DPEP/curled_lib/script/input_params.csv), and output pdb file is named as `H_<len>_<phi0>_<delphi>_<psi0>_<delpsi>_<phase>.pdb`.
We also provide helix scaffold library at various lengths (28 aa and 35 aa) already generated in this work. Click [here](https://pan.baidu.com/s/1lKv6-XoMh6dJfG7dW_JR4A?pwd=14kv) (extraction code: 14kv) to download.

### Tutorial for helix scaffold fitting
#### Flip the target segment into D-type
Before docking, please flip your target segment to D-type, with residue names unchanged. Note that input file type should be a pdb file **with hydrogens removed**.
```shell
chmod +x $DPEP/docking/mirror_target/mirror_target.sh
$DPEP/docking/mirror_target/mirror_target.sh -i your_input_file.pdb -o your_output_file.pdb
```
The default output of `-o` is your_input_file_mirror.pdb

#### L- scaffold fitting to D- reference segment
You can fit scaffolds to extracted D- reference helical segment.
```shell
gcc $DPEP/docking/ScaffoldFitting/FitequationD.c -o $DPEP/docking/ScaffoldFitting/FitequationD -lm
```
Running FitequationD to generate fitted scaffolds:
```shell
$DPEP/docking/ScaffoldFitting/FitequationD [Dhelixtemplate.pdb] [fitted_output.pdb] [startResidueForFit default:0]
```
where `Dhelixtemplate.pdb` refers to the flipped target segment, `fitted_output.pdb` is the file name for fitted results. 'startResidueForFit' defines the start residue for fitting.





### Tutorial for HelixScaffoldDocking
#### Flip the target into D-type
Before docking, please flip your target to D-type as mentioned above.
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

#### Docking process
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

The input scaffold file and output file are separated by spaces.


## About database
We provide information of helix ligand-protein complex database in `$DPEP/database/Hlig_protein_database.csv`, which contains PDBID, residue numbers of helical ligands and fitting results.

## License

This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.

## Acknowledgements
This design procedure involves multiple softwares related to protein design. We acknowledge and thank the developers of Naccess, Rosetta, gromacs, gmx_MMPBSA, AggreScan and PASTA2.0 for their incredible and hard work.

