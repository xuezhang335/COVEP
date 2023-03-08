# COVEP -- Pepline for B/T cell epitopes prediction

# Introduction
COVEP is an ensemble pipeline for coronavirus B/T cell epitopes prediction, which provides [web sever](https://pgx.zju.edu.cn/coviedb2/tools/COVEP) , [Docker image](https://hub.docker.com/repository/docker/biopharm/COVEP/general) and also supports downloading from [Github-COVEP](https://github.com/zjupgx/COVEP). It has been used for SARS-CoV-2, SARS-CoV, MERS-CoV and other four human coronaviruses (hCoV-229E, hCoV-HKU1, hCoV-OC43, hCoV-NL63) epitopes prediction and shared conserved epitopes analysis for peptide-based vaccine design and the result are presented in our database [COVIEdb](http://pgx.zju.edu.cn/coviedb/).

Contact: zhanzhou@zju.edu.cn
<br><br>

# Download and installation
COVEP is a python program, and requires python to be installed on your system. We highly recommend using [Conda](https://www.anaconda.com/products/distribution) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) to manage environments, i.e. use Conda to create an environment with dependencies installed where you can use COVEP. It’s possible to use other methods, but this documentation presupposes that you have Conda installed.

It is worth noting that if you use docker to run, then you do not need to consider the environment, all dependencies have been packaged in the docker image, see [section-Docker](#docker) for instructions.

## 1. Git
COVEP runs under Linux System, and all the dependencies should be properly pre-installed. The file of `cov.yaml` lsits the version of the dependency packages.

Download the latest version from [GitHub-COVEP](https://github.com/zjupgx/COVEP) :

    git clone git@github.com:zjupgx/COVEP.git

Unzip the source code and go into the directory:

    tar xvzf cov-*.tar.gz
    cd cov

Create conda environment:

    conda create -f cov.yaml

Modify the arguemrnts of `seqfile`, `hlafile`, `outdir` and invoke the shell script. See [General Usage: parameters](#general-usage) for parameter description:

    bash cov.sh

    ```
    #! /bin/bash

    seqfile=/root/cov/test/data/2019-ncov-test.txt
    hlafile=/root/cov/test/data/hla_allele-test.csv
    outdir=/root/cov/test/output

    python /root/cov/scripts/main.py -bl all -t1 mhcflurry,netmhcpan -t2 all -s $seqfile -hlafile $hlafile -o $outdir
    ```
You can also directly run the `main.py` python file in terminal.


## 2. Docker
This approach requires [Docker](https://docs.docker.com/) to be installed on your system. Conda is installed in the image, so it will be relatively large, about 6.3G. For easy editing and visualization, jupyter has been built into the image runtime environment. The workdir of COVEP is `/root/cov`.

Pull the [cov](https://hub.docker.com/repository/docker/zhangxue123/cov) image from dockerhub:

    docker pull biopharm/cov:v1.0

run image and enter the generated container in bash mode:

    docker run -it -p 8888:8888 --name [your/container/name] biopharm/cov:v1.0 bash

you can start a Jupyter Notebook server and interact via your browser:

    conda activate cov
    jupyter notebook --notebook-dir=/root/cov --ip='*' --no-browser --allow-root

Testing in the container, and modify the `seqfile`, `hlafile`, `outdir` in the cov.sh:

    bash cov.sh
    
    ```
    #! /bin/bash

    seqfile=/root/cov/test/data/2019-ncov-test.txt
    hlafile=/root/cov/test/data/hla_allele-test.csv
    outdir=/root/cov/test/output

    python /root/cov/scripts/main.py -bl all -t1 mhcflurry,netmhcpan -t2 all -s $seqfile -hlafile $hlafile -o $outdir
    ```

## 3. Web sever
COVEP provides a [web sever](https://pgx.zju.edu.cn/coviedb2/tools/COVEP) for the convenience of users. It only needs to submit sequences and HLA genes to make predictions without configuring the environment.
<br><br>

# General Usage

COVEP needs to be called by running the `main.py` python file with inputting the specified parameters. Users can modify the parameters of the `bash.cov` script to run, or run directly on the terminal command line.

***Running command***

    python main.py -bl iedb7 -t1 netmhcpan,mhcflurry -t2 all -s path/to/your/input/seqfile -hlafile path/to/your/hlafile -o directory/to/store/predicted/results

***parameters***

    -bl         Tools used for linear B cell epitopes prediction, use `iedb7` to use seven IEDB linear B cell epitope tools.
    -t1         Tools used for MHC class I T cell epitopes prediction, choose from [netmhcpan,mhcflurry,deephlapan], eg. netmhcpan,mhcflurry. Use `all` to select all available tools.
    -t2         Tools used for MHC class II T cell epitopes prediction, choose from [netmhc2pan,mixmhc2pred], eg. netmhc2pan. Use `all` to select all available tools.
    -seqfile    txt file contains one or more fasta sequences, see test/data/2019-ncov-test.txt for example.
    -hla        str, this is necessary if t1 or t2 is selected, such as `HLA-A11:01,HLA-A24:02,` for t1 or `DRB1_01_01,DRB1_01_02` for t2'.
    -hlafile    csv file, `-hlastr` will be covered if this option is specified. This is necessary if t1 and t2 both are selected, see test/data/hla_allele.csv for example.
    -o          Directory to store predicted results. User must have write privilege. If omitted, the current directory will be applied.

***Supported prediction tools***

    | Epitope | Length | Tools | 
    | :----- | :----- | :----- |
    | Linear B | pisition score | Chou-Fasman, Emini, Karplus-Schulz, <br />Kolaskar-Tongaonkar, Parker, Bepipred |
    | MHC class I  | 8-14mer | DeepHLApan, NetMHCpan-4.1, MHCflurry |
    | MHC class II | 12-25mer | NetMHCIIpan-4.0, MixMHC2pred |

Supported prediction tools are shown as in the table. Please note that linear B cell epitopes prediction is depend on the network, and DeepHLApan and NetMHC2pan maybe slow.

***HLA allele and threshold***

Supported allele and threshold of prediction tools seen in the file `HLA_Support.xlsx`
<br><br>
# Update log

### 2022.04
***v1.0***

1. Support for linear B cell epitopes, class I and class II T cell epitopes prediction.

### 2022.10
***v1.1***
1. For t1 and t2 prediction, predicted score is added to results.
2. Support the use of single or several predicting tools. 


