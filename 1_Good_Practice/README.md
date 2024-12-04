# Good Practice

## Connect to 1 HPC/cluster to handle the big dataset (optional) 
Each experiment will acquire different temporaly memories of the computer. I decide to use **plaFRIM** to handle my experiment but you can use other cluster/HPC or eventually your local computer.
Ensure to upload your folders in the cluster

```
$ scp -r /path/to/your/file/citeseq_wk4 /path/to/your/cluster/repository/
```
Then sign in with ssh

```
$ ssh myusername@plafrim
```
Navigate in the your folder

```
[my_username@plafrim]$ cd /path/to/citeseq_wk4
```

## Create an environment dedicated to downstream analysis
An specific environment is a good way to handle dependencies of the packages/tools that you will use
I used conda 

