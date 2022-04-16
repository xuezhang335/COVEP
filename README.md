# Pepline for B/T cell epitopes prediction 

## Despriction

This is an ensemble pipeline for coronavirus B/T cell epitopes prediction, which provides a web sever and also supports downloading to run locally. It has been used for SARS-CoV-2, SARS-CoV, MERS-CoV and other four human coronaviruses (hCoV-229E, hCoV-HKU1, hCoV-OC43, hCoV-NL63) epitopes prediction and shared conserved epitopes analysis for peptide-based vaccine design and the result are presented in our database.  [COVIEdb – A database for potential immune epitopes of coronaviruses](http://biopharm.zju.edu.cn/coviedb/)

Workflow.

Contact: zhanzhou@zju.edu.cn

## 	Web Sever

[COVIEdb – A database for potential immune epitopes of coronaviruses](http://biopharm.zju.edu.cn/coviedb/)



## Download and installation

### 1. Docker (Recommend)

The Installation of Docker can be seen in https://docs.docker.com/

Pull the image of cov from dockerhub:

```sh
docker pull biopharm/cov:v1.0
```

run the image in bash mode:

```sh
 docker run -it -p 8888:8888 --name [your/container/name] biopharm/cov:v1.0 bash
```

you can start a Jupyter Notebook server and interact  via your browser:

```sh
conda activate cov
jupyter notebook --notecook-dir=/root/cov --ip='*' --no-browser --allow-root
```

test in the container, and modify the seqfile, hlafile, outdir in the cov.sh:

```sh
bash cov.sh

'''
#! /bin/bash

seqfile=/root/cov/test/data/2019-ncov-test.txt
hlafile=/root/cov/test/data/hla_allele-test.csv
outdir=/root/cov/test/output

python /root/cov/scripts/main.py -bl -t1 -t2 -s $seqfile -hlafile $hlafile -o $outdir
'''

```



### 2. Git (All the dependencies should be properly installed)

#### System

Linux

#### Dependencies

See cov.yaml

#### Steps

Download the latest version of cov from https://github.com/xuezhang335/CoVEP
    

```shell
git clone https://github.com/jiujiezz/cov.git
```

Unzip the source code and go into the directory by using the following command:

```shell
tar xvzf cov-*.tar.gz

cd cov
```

Invoke the shell script:

```shell
bash cov.sh
```



## General Usage



```shell
Usage    
python main.py --bl --t1 --t2 -s path/to/your/input/seqfile -hlafile path/to/your/hlafile -o directory/to/store/predicted/results

optional arguments:
--bl     Prediction for linear B cellepitopes.           store_true
--t1     Prediction for MHC class I T cell  epitopes.       store_true
--t2     Prediction for MHC class II T cell  epitopes.      store_true
-seqfile   txt file contains one or more fasta sequences.
-hlafile   csv file, which is necessary if t1 and t2 both are selected.
-hla str,  is necessary if t1 or t2 is selected.
-o   Directory to store predicted results.


Supported length
B epitopes：pisition score；
T epitopes class I：8-11mer；
T epitopes class II：15-25mer；

Supported allele and threshold of prediction tools seen in the file (HLA_Support.xlsx)
```



| Epitope      | Tools                                                        |
| :----------- | ------------------------------------------------------------ |
| Linear B     | Chou-Fasman, Emini, Karplus-Schulz,  <br />Kolaskar-Tongaonkar, Parker, Bepipred |
| MHC class I  | DeepHLApan, NetMHCpan-4.1, MHCflurry                         |
| MHC class II | NetMHCIIpan-4.0, MixMHC2pred                                 |

Prediction tools are shown as in the table, please note that linear B cell epitopes prediction is depend on the network and NetMHCIIpan maybe slow when meets long sequence.





## Update log

### 2022.04

v1.0
Support for linear B cell epitopes, class I and class II T cell epitopes prediction.



