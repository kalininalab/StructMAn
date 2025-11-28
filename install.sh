#!/bin/bash

# Absolute path to this script.
SCRIPT=$(readlink -f $0)
# Absolute path this script is in.
SCRIPTPATH=$(dirname "$SCRIPT")

all_good=true

# Stop script on errors
#set -euo pipefail

#Init constants
current_python_version="3.12"
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
path_to_dssp_build_script="$SCRIPTPATH"/build_smssp/build.sh
path_to_dssp_binary="$SCRIPTPATH"/smssp

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

verbose_stdout=3
verbose_stderr=4

if [ "$verbose" = true ]; then
    eval "exec $verbose_stdout>&1"
    eval "exec $verbose_stderr>&2"
else
    eval "exec $verbose_stdout>/dev/null"
    eval "exec $verbose_stderr>/dev/null"
fi

#Check if conda environment already exits, create it if not
env_list_result=$(conda env list | grep -w "$env_name")
if [ -z "$env_list_result" ]
then
    echo "Conda environment with name $env_name not in current environment list, setting up new environment ..."
    if ! conda create -n "$env_name" python=$current_python_version -y >&$verbose_stdout; then
        echo "Installation failed"
        exit 1
    fi
else
    echo "Conda environment with name $env_name already in environment list."
fi

#Locate conda source to activate conda inside bash script shell
conda_base_path=$(conda info --base)
conda_bash_path="$conda_base_path"/etc/profile.d/conda.sh

current_env=$(conda info | grep 'active environment' | awk '{sub(/^[ \t\r\n]+active environment.:./,""); print}')
if ! [ "$current_env" = "base" ]
then
    echo "Conda needs to be in base environment, but is $current_env"
    echo "Please call 'conda deactivate'."
    echo "exiting installer ..."
    exit 1
fi


#activate the environment
{
    echo "Activating environment inside shell ..."
    source "$conda_bash_path"

    new_env_path=$(conda env list | awk -v name="$env_name" '/^[^#]/{ if ($1 == name) {print $2} }')

    #Installing custom DSSP binary

    if ! $path_to_dssp_build_script; then
        echo "Installation failed in building smssp"
        exit 1
    fi
    if ! mv "$path_to_dssp_binary" "$new_env_path"/bin/smssp; then
        echo "Installation failed"
        exit 1
    fi

    conda activate "$env_name"
    echo "$env_name" activated
} >&$verbose_stdout

mamba_version_test_output=$(mamba --version 2>/dev/null)

if [ -z "$mamba_version_test_output" ]
then
    {
        conda install -y -c conda-forge mamba 
        echo "Mamba installed"
    } >&$verbose_stdout 2>&$verbose_stderr
else
    echo "Mamba already installed: $mamba_version_test_output" >&$verbose_stdout
fi

mamba update pip -y -c conda-forge >&$verbose_stdout

#install dependencies
{
    echo "Installing package cairo ..."
    mamba install -y cairo==1.18.4 -c conda-forge
    export PKG_CONFIG_PATH="$new_env_path"/envs/"$env_name"/lib/pkgconfig:"$PKG_CONFIG_PATH"
    echo "Installing package pycairo ..."
    if ! mamba install -y pycairo==1.28.0 -c conda-forge; then
        echo "Installation failed"
        exit 1
    fi
    echo "Installing package importlib-metadata ..."
    mamba install -y importlib-metadata
    echo "Installing package keyring ..."
    if ! mamba install -y keyring -c conda-forge; then
        echo "Installation failed"
        exit 1
    fi
    echo "Installing package pkginfo ..."
    mamba install -y pkginfo
    echo "Installing package requests-toolbelt ..."
    mamba install -y requests-toolbelt
    echo "Installing package pymol ..."
    mamba install -y -c conda-forge pymol-open-source
    echo "Installing package gawk ..."    
    mamba install -y gawk
    echo "Installing package sqlite ..."
    if ! mamba install -y -c conda-forge libsqlite --force-reinstall; then
        echo "Mamba failed, try to fallback to conda..."
        if ! conda install -y -c conda-forge libsqlite --force-reinstall; then
            echo "Installation failed"
            exit 1
        fi
    fi
    mamba update -y -c conda-forge libsqlite
    if ! mamba install -y sqlite -c conda-forge --force-reinstall; then
        echo "Installation failed"
        exit 1
    fi
    mamba update -y sqlite -c conda-forge
    echo "Installing package openbabel ..."    
    mamba install -y openbabel -c conda-forge
} >&$verbose_stdout


spectrumany_target_path=$"$SCRIPTPATH/structman/scripts/spectrumany.py"
if [ ! -f "$spectrumany_target_path" ]
then
    {
        pushd "$SCRIPTPATH/structman/scripts/"
        wget https://raw.github.com/Pymol-Scripts/Pymol-script-repo/master/spectrumany.py
        popd
    }
fi


#install the main package
echo "Installing StructMAn source code using pip ..."
if ! pip install "$SCRIPTPATH" >&$verbose_stdout; then
    echo "Installation failed"
    exit 1
fi

#install mmseqs2
echo "Installing MMseqs2 ..."
mamba install -y -c bioconda -c conda-forge mmseqs2==15.6f452 >&$verbose_stdout

#install wkhtmltopdf
mamba install -y -c conda-forge wkhtmltopdf >&$verbose_stdout

if [ -z $storage_folder ]
then
    storage_folder="$new_env_path"/share/structman
    mkdir "$storage_folder"
fi

if [ -d $storage_folder ]
then
    storage_folder=$(realpath "$storage_folder")
else
    mkdir "$storage_folder"
    storage_folder=$(realpath "$storage_folder")
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

if [ -z $local_database_folder ]
then
    local_database_folder="$storage_folder"
fi

if [ -d "$local_database_folder" ]
then
    local_database_folder=$(realpath "$local_database_folder")
else
    mkdir "$local_database_folder"
    local_database_folder=$(realpath "$local_database_folder")
fi

#Installing specific boost
if ! mamba install -y zlib==1.3.1 -c conda-forge >&$verbose_stdout; then
    echo "Installation failed"
    exit 1
fi

if ! mamba install -y libfreetype6==2.13.3 -c conda-forge >&$verbose_stdout; then
    echo "Installation failed"
    exit 1
fi

if ! mamba install -y libboost==1.84 -c conda-forge >&$verbose_stdout; then
    echo "Installation failed"
    exit 1
fi



#Installing grpc
pip install grpcio >$verbose_stdout

#Init config file
structman_config_path="$new_env_path"/lib/python"$current_python_version"/site-packages/structman/structman_config.txt
cp "$SCRIPTPATH/config_template.txt" "$structman_config_path"
echo "Setting up personal structman config file:"
echo "    $structman_config_path"
{
    if ! structman config mmseqs_tmp_folder "$tmp_folder_path" -c "$structman_config_path"; then
        echo "Installation failed"
        exit 1
    fi
    structman config dssp_path smssp -c "$structman_config_path"
    structman config mmseqs2_db_path "$storage_folder"/pdbba_mmseqs2_search_db
} >&$verbose_stdout 2>&$verbose_stderr

#install modeller
if ! [ -z "$modeller_key" ]
then
    echo "Installing Modeller ..."
    {
        echo "Installing modeller package ..."
        mamba install -y modeller -c salilab
        echo "Setting the given modeller key to the modeller config ..."
        structman set_m_key "$modeller_key"
    } >&$verbose_stdout 2>&$verbose_stderr
fi

echo "Setting StructMAn config and database ..."

if ! [ -z "$database_account_name" ]
then
    database_name="$database_account_name"_main_db
    {
        echo "Configuring the database connection ..."
        structman config db_user_name "$database_account_name" -c "$structman_config_path"
        structman config db_password "$database_password" -c "$structman_config_path"
        structman config db_name "$database_name" -c "$structman_config_path"
        structman config db_address "$database_server" -c "$structman_config_path"
        echo "Setting up the database ..."
        structman database create --compute_ppi
    } >&$verbose_stdout 2>&$verbose_stderr
else
    database_name=main_db
    local_database_path="$local_database_folder"/local_structman_database.db
    {
        echo Setting up local database here: "$local_database_path"
        if ! structman config local_db_path "$local_database_path" -c "$structman_config_path"; then
            echo "Installation failed"
            exit 1
        fi
        structman config db_name "$database_name" -c "$structman_config_path"
        structman config db_address "-"
        structman database create --compute_ppi
    } >&$verbose_stdout 2>&$verbose_stderr
fi

structman update check_search_db >&$verbose_stdout 2>&$verbose_stderr

if $all_good;
then
    echo "StructMAn successfully installed, please activate the right conda environment before using it:"
    echo "    conda activate $env_name"
    exit 0
else
    echo "Installation failed"
    exit 1
fi