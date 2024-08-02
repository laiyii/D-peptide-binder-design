# Targeting SARS-CoV-2 Main Protease with D-peptides
Welcome to D-peptide binder design project! This repository is the source code for Targeting SARS-CoV-2 Main Protease with D-peptides.
![workflow](https://github.com/laiyii/D-peptide-binder-design/blob/main/Dpep_fig1.png)
If you have any questions, please discuss in [Issues](https://github.com/laiyii/D-peptide-binder-design/issues).

## Installation




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
We also provide helix scaffold library with various lengths (21 aa, 24 aa, 28 aa, 31 aa, 35 aa, 38 aa, and 42 aa) already generated in this work:
```shell
$HSD/curled_lib/generated_lib
```

### Docking of the helical scaffolds to the target
#### Flip the target into D-type
Before docking, please flip your target to D-type. Note that input file type should be pdb **with hydrogens removed**.
```shell
chmod +x $HSD/docking/mirror_target/mirror_target.sh
$HSD/docking/mirror_target/mirror_target.sh -i your_input_file.pdb -o your_output_file.pdb
```
The default output of `-o` is your_input_file_mirror.pdb

#### Helix scaffolds docking


### Loop modeling with CCD

### Sequence design


