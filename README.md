# StructMAn

**Structural Mutation Annotation**

StructMAn is a computational pipeline and tool kit to annotate protein sequences to the protein structure level combined with an automatical comprehensive structural analysis. The method was initially designed to focus on specifically given individual positions in protein sequences, but evolved over time to process full protein sequences.  \
Apart from the short overview, you can find a more detailed explanation of the tool on 
[ReadTheDocs](https://structman.readthedocs.io/en/main/). 

## Installation

StructMAn needs to be deployed in a a conda environment, therefore a installation of conda (or miniconda or mamba) is a dependency.
>[!TIP]
>For example, miniconda can be installed following this [website](https://docs.anaconda.com/miniconda/miniconda-install/)

* Download structman repo
  ```
  git clone https://github.com/kalininalab/StructMAn.git
  ```
* Navigate into repo folder
  Cloning the repo will create a folder named `StructMAn` that contains the installation script `install.sh`
  ```
  cd StructMAn
  ```

* Install structman (installation script inside the cloned repo)


>[!WARNING]
>**You need to be in conda base environment to call the install.sh script.**
  
  ```
  ./install.sh [optional flags listed below]
  ```
  * `-e 'environment name'` Name of the conda environment StructMAn should installed into. Default: 'structman'
  * `-d 'database account name'` Account name to access our MySQL database services. Should always provided together with a password (`-p 'database password'`) and MySQL database server connection (`-c 'database_server'`). **If no database server available, use the local sqlite3 version (`-f`)**
  * `-p 'database password'` Password (corresponded to the account name) to access the MySQL database services.
  * `-c 'database_server'` MySQL database server connection.
  * `-m 'Modeller license key'` Modeller license key phrase. **Only required when StructMAn is used for automated Modelling.** If not given, modeller will not be installed and the modelling functionalities of StructMAn will not be usable. If a wrong key is given, modeller will be installed, but it will fail on usage. (Granted for academic use only by [salilab](https://salilab.org/modeller/))
  * `-s 'storage folder'` Path to a folder, where structman may write files and store some necessary data. On installation a mmseqs2 indices table with around 1.5Gb will be stored there. On runtime, temporary files will be written there. Default: /share/structman/ in the corresponding environment folder
  * `-f 'local database folder'` Path to a folder, where structman stores a sqlite3 database. Not necessary when a database connection via `-c`, `-d`, and `-p` were given. The size of the database depends on the size of processed inputs and can reach several hundred Gb. Default: the storage folder given with `-s`
  * `-v` Adding verbose Output. **When having issues with StructMAn installation** please provide verbose Output. It is quite large, so best redirect it with `&> structman_installer_error.txt`

 The installation script will generate a structman config file that will be used automatically. This file can be copied/modified for further usage (for example using specific database instances).

 * Switch into the corresponding environment
   ```
   conda activate 'environment name'
   ```
 * Run StructMAn
   ```
   structman -i my_example.fasta
   ```
