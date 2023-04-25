# HSHMP_2022

Code for reproducing the figures and results in the preprint "Accurate quantification of single-nucleus and single-cell RNA-seq transcripts" by Kristján Eldjárn Hjörleifsson, Delaney K. Sullivan, Guillaume Holley, Páll Melsted and Lior Pachter

## Introduction

Please follow the steps below in order to reproduce the results of the preprint. Set all the paths to be relative to the directory HSHMP_2022.

<pre>main_path="$(pwd)/HSHMP_2022"
kallisto="$main_path/kallisto-D/build/src/kallisto"
bustools="$main_path/bustools/build/src/bustools"
cellranger3="$main_path/cellranger/cellranger-3.1.0/cellranger"
cellranger4="$main_path/cellranger/cellranger-4.0.0/cellranger"
cellranger5="$main_path/cellranger/cellranger-5.0.1/cellranger"
cellranger6="$main_path/cellranger/cellranger-6.1.2/cellranger"
cellranger7="$main_path/cellranger/cellranger-7.0.1/cellranger"
salmon="$main_path/salmon-latest_linux_x86_64/bin/salmon"
</pre>

## Download software

### Kallisto-D

commit 0c80683b37b651c93d67fc8cf02658d4ab2cd2a3

<pre>cd $main_path
rm -rf kallisto-D
git clone https://github.com/pachterlab/kallisto-D
cd kallisto-D && mkdir -p build && cd build
cmake .. -DZLIBNG=ON && make</pre>

### Bustools

commit f4fd12a5205772eb7e62a04a7ebd8b40835805b8

<pre>cd $main_path
rm -rf bustools
git clone -b dlist https://github.com/BUStools/bustools
cd bustools && mkdir -p build && cd build
cmake .. && make</pre>

### kb-python

commit b0700bb78625bc82918396d86bcd538f8d9c838c

<pre>cd $main_path
yes|pip uninstall kb-python
pip install --user git+https://github.com/pachterlab/kb_python@dlist</pre>

### CellRanger

Note: CellRanger needs to be installed manually. Versions are as follows:

* Cell Ranger v3.1.0 (Released July 24, 2019. Downloaded November 22, 2022)
* Cell Ranger v4.0.0 (Released July 7, 2020. Downloaded November 22, 2022)
* Cell Ranger v5.0.1 (Released December 16, 2020. Downloaded November 22, 2022)
* Cell Ranger v6.1.6 (Released October 25, 2021. Downloaded October 7, 2022)
* Cell Ranger v7.0.1 (Released August 18, 2022. Downloaded October 7, 2022)

### salmon-alevin-fry

salmon version 1.10.0; alevin-fry version 0.8.0; pyroe 0.6.4

<pre>cd $main_path
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && tar -xzvf salmon-1.10.0_linux_x86_64.tar.gz
export RUSTUP_HOME=./.rustup/
export CARGO_HOME=./.cargo/
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
./.cargo/bin/cargo install alevin-fry
pip install pyroe</pre>


## STARSolo simulations

Navigate to sim_starsolo and run the scripts there

