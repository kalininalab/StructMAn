# StructMAn

**Structural Mutation Annotation**

StructMAn is a computational pipeline and tool kit to annotate protein sequences to the protein structure level combined with an automatical comprehensive structural analysis. The term "Mutation" is rather historical due the method was initially designed to focus on specifically given individual positions in protein sequences, but evolved over time to process full protein sequences.  \
Apart from the short overview, you can find a more detailed explanation of the tool on 
[ReadTheDocs](https://structman.readthedocs.io/en/main/). 

## Installation

StructMAn needs to be deployed in a a conda environment, therefore a installation of conda (or miniconda or mamba) is a dependency.

* Download structman repo
  ```
  git clone https://github.com/kalininalab/StructMAn.git
  ```

* Install structman (installation script inside the cloned repo)
  ```
  ./install.sh [-e 'environment name'] [-d 'database account name'] [-p 'database password'] [-c 'database_server'] [-m 'Modeller license key'] [-s 'storage folder'] [-f 'local database folder']
  ```
  * `-e` Name of the conda environment StructMAn should installed into. Default: 'structman'
  * `-d` Account name to access our MySQL database services. Should always provided together with a password (`-p`) and MySQL databaser server connection (`-c`).
  * `-p` Password (corresponded to the account name) to access the MySQL database services.
  * `-m` Modeller license key phrase. If not given, modeller will not be installed and the modelling functionalities of StructMAn will not be usable. If a wrong key is given, modeller will be installed, but it will fail on usage. (Granted for academic use only by [salilab](https://salilab.org/modeller/))
  * `-s` Path to a folder, where structman may write files and store some necessary data. On installation a mmseqs2 indices table with around 1.5Gb will be stored there. On runtime, temporary files will be written there. Default: /share/structman/ in the corresponding environment folder
  * `-f` Path to a folder, where structman stores a sqlite3 database. Not necessary when a database connection via `-c`, `-d`, and `-p` were given. The size of the database depends on the size of processed inputs and can reach several hundred Gb. Default: the storage folder given with `-s`

 The installation script will generate a structman config file that will be used automatically. This file can be copied/modified for further usage (for example using specific database instances).

 * Switch into the corresponding environment
   ```
   conda activate 'environment name'
   ```
 * Run StructMAn
   ```
   structman -i my_example.fasta
   ```
