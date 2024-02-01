#!/bin/bash

# Absolute path to this script.
SCRIPT=$(readlink -f $0)
# Absolute path this script is in.
SCRIPTPATH=$(dirname "$SCRIPT")

#Init constants
current_python_version="3.11"
username=$(whoami)

#Init default arguments
env_name=structman
modeller_key=""
database_account_name=""
database_password=""
database_server=""
local_database_folder=""
storage_folder=""
verbose=false

installer_temp_folder=$(mktemp -d -t structman.XXXXXXX)
trap 'rm -rf -- "$installer_temp_folder"' EXIT
path_to_dssp_patch="$SCRIPTPATH"/structman/scripts/dssp.patch

#Parse arguments
while getopts e:m:d:p:f:s:c:v flag
do
    case "${flag}" in
        e) env_name=${OPTARG};;
        m) modeller_key=${OPTARG};;
        d) database_account_name=${OPTARG};;
        p) database_password=${OPTARG};;
        f) local_database_folder=${OPTARG};;
        s) storage_folder=${OPTARG};;
        c) database_server=${OPTARG};;
        v) verbose=true;;
    esac
done

if [ -z "$database_account_name" ] && [ -z "$local_database_folder" ]
then
    echo "ERROR - Invalid Input"
    echo "You need to provide either a database connection information or a path to a folder to store a local database"
    echo "Please provide a local database folder with -f [path to folder] or"
    echo "a database account name with -d [account name] and -p [password]"
    exit 1
fi

#Check if conda environment already exits, create it if not
env_list_result=$(conda env list | grep "$env_name")
if [ -z "$env_list_result" ]
then
    echo "Conda environment with name $env_name not in current environment list, setting up new environment ..."
    if $verbose
    then
        conda create -n "$env_name" python=$current_python_version -y
    else
        conda create -n "$env_name" python=$current_python_version -y >/dev/null
    fi
else
    echo "Conda environment with name $env_name already in environment list."
fi

#Locate conda source to activate conda inside bash script shell
conda_base_path=$(conda env list | grep base | sed -e 's/.* \([^ ]*\)$/\1/g') # Regex fetches the base path to conda
conda_bash_path="$conda_base_path"/etc/profile.d/conda.sh

new_env_path=$(conda env list | grep "$env_name" | sed -e 's/.* \([^ ]*\)$/\1/g')

#activate the environment
if $verbose
then
    echo "Activating environment inside shell ..."
    source "$conda_bash_path"
    conda activate "$env_name"
    echo "$new_env_path"" activated"
else
    source "$conda_bash_path" >/dev/null
    conda activate "$env_name" >/dev/null
fi

mamba_version_test_output=$(mamba --version 2>/dev/null)

if [ -z "$mamba_version_test_output" ]
then
    if $verbose
    then
        conda install -y -c conda-forge mamba
    else
        conda install -y -c conda-forge mamba &> /dev/null
    fi
fi

#install dependencies
if $verbose
then
    echo "Installing package cairo ..."
    mamba install -y cairo
    export PKG_CONFIG_PATH="$new_env_path"/envs/"$env_name"/lib/pkgconfig:"$PKG_CONFIG_PATH"
    echo "Installing package pycairo ..."
    mamba install -y pycairo    
    echo "Installing package importlib-metadata ..."
    mamba install -y importlib-metadata
    echo "Installing package keyring ..."
    mamba install -y keyring
    echo "Installing package pkginfo ..."
    mamba install -y pkginfo
    echo "Installing package requests-toolbelt ..."
    mamba install -y requests-toolbelt
    echo "Installing package pymol ..."    
    mamba install -y -c conda-forge pymol-open-source
else
    mamba install -y cairo >/dev/null
    export PKG_CONFIG_PATH="$new_env_path"/envs/"$env_name"/lib/pkgconfig:"$PKG_CONFIG_PATH"
    mamba install -y pycairo >/dev/null
    mamba install -y importlib-metadata >/dev/null
    mamba install -y keyring >/dev/null
    mamba install -y pkginfo >/dev/null
    mamba install -y requests-toolbelt >/dev/null
    mamba install -y -c conda-forge pymol-open-source >/dev/null
fi

#install the main package
echo "Installing StructMAn source code using pip ..."
if $verbose
then
    pip install "$SCRIPTPATH"
else
    pip install "$SCRIPTPATH" >/dev/null
fi

#install mmseqs2
echo "Installing MMseqs2 ..."
if $verbose
then
    mamba install -y -c bioconda mmseqs2
else
    mamba install -y -c bioconda mmseqs2 >/dev/null
fi

if [ -z $storage_folder ]
then
    storage_folder="$new_env_path"/share/structman
    mkdir "$storage_folder"
fi

if ! [ -d "$storage_folder/$username" ]
then
    mkdir "$storage_folder/$username"
    echo "Needed to create a user-specific folder on scratch to locate the tmp folder:"
    echo "   $storage_folder/$username"
fi

tmp_folder_path="$storage_folder/$username/tmp"
if ! [ -d "$tmp_folder_path" ]
then
    mkdir "$tmp_folder_path"
    echo "Needed to create the tmp folder:"
    echo "   $tmp_folder_path"
fi

#Fetching, patching, and installing custom DSSP binary
pushd "$installer_temp_folder"
wget https://github.com/cmbi/hssp/releases/download/2.2.8/xssp-2.2.8.tar.gz 
tar -zxf xssp-2.2.8.tar.gz
cd xssp-2.2.8
patch -p 1 < "$path_to_dssp_patch"
./configure
make mkdssp
mv mkdssp "$new_env_path"/bin/smssp
popd

#Init config file
structman_config_path="$new_env_path"/lib/python"$current_python_version"/site-packages/structman/structman_config.txt
cp "$SCRIPTPATH/config_template.txt" "$structman_config_path"
echo "Setting up personal structman config file:"
echo "    $structman_config_path"
if $verbose
then 
    structman config mmseqs_tmp_folder "$tmp_folder_path" -c "$structman_config_path"
    structman config dssp_path smssp -c "$structman_config_path"
    structman config mmseqs2_db_path "$storage_folder"/pdbba_search_db_mmseqs2
else
    structman config mmseqs_tmp_folder "$tmp_folder_path" -c "$structman_config_path" &>/dev/null
    structman config dssp_path smssp -c "$structman_config_path" &>/dev/null
    structman config mmseqs2_db_path "$storage_folder"/pdbba_search_db_mmseqs2 &>/dev/null
fi

#install modeller
if ! [ -z "$modeller_key" ]
then
    echo "Installing Modeller ..."
    if $verbose
    then
        echo "Adding salilab to the conda channels ..."
        mamba config --add channels salilab
        echo "Installing modeller package ..."
        mamba install -y modeller
        echo "Setting the given modeller key to the modeller config ..."
        structman set_m_key "$modeller_key"
    else
        mamba config --add channels salilab &>/dev/null
        mamba install -y modeller >/dev/null
        structman set_m_key "$modeller_key" &>/dev/null
    fi
fi

echo "Setting StructMAn config and database ..."

if ! [ -z "$database_account_name" ]
then
    database_name="$database_account_name"_main_db
    if $verbose
    then
        echo "Configuring the database connection ..."
        structman config db_user_name "$database_account_name" -c "$structman_config_path"
        structman config db_password "$database_password" -c "$structman_config_path"
        structman config db_name "$database_name" -c "$structman_config_path"
        structman config db_address "$database_server" -c "$structman_config_path"
        echo "Setting up the database ..."
        structman database create --compute_ppi
    else
        structman config db_user_name "$database_account_name" -c "$structman_config_path" &>/dev/null
        structman config db_password "$database_password" -c "$structman_config_path" &>/dev/null
        structman config db_name "$database_name" -c "$structman_config_path" &>/dev/null
        structman config db_address "$database_server" -c "$structman_config_path" &>/dev/null
        structman database create --compute_ppi &>/dev/null
    fi
else
    database_name=main_db
    local_database_path="$local_database_folder"/local_structman_database.db
    if $verbose
    then
        echo Setting up local database here: "$local_database_path"
        structman config local_db_path "$local_database_path" -c "$structman_config_path"
        structman config db_name "$database_name" -c "$structman_config_path"
        structman config db_address "-"
        structman database create --compute_ppi
    else
        structman config local_db_path "$local_database_path" -c "$structman_config_path" &>/dev/null
        structman config db_name "$database_name" -c "$structman_config_path" &>/dev/null
        structman config db_address "-" &>/dev/null
        structman database create --compute_ppi &>/dev/null
    fi
fi

if $verbose
then
    structman update check_search_db
else
    structman update check_search_db &>/dev/null
fi

echo "StructMAn successfully installed, please activate the right conda environment before using it:"
echo "    conda activate $env_name"

exit 0