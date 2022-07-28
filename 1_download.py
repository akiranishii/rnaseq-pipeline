#!/bin/python

'''
Downloading data from Beijing Genomics Institute
'''

import re
import subprocess
import os

def download(file):
    '''Download a file using wget'''
    if file.count("/") > 0:
        subdir = re.findall(".+(?=/)", file)[0]
        call = f"wget -P {subdir} {path}/{file}"
    else:
        call = f"wget {path}/{file}"
    print(call)
    subprocess.call(call, shell=True)

if __name__ == '__main__':
    os.mkdir("F21FTSUSAT0328_MOUhnjzR")
    os.chdir("F21FTSUSAT0328_MOUhnjzR")
    path = "ftp://20210429F21FTSUSAT0328:MOUhnjzRlqp_0429@cdts-hk.genomics.cn/F21FTSUSAT0328_MOUhnjzR"
    # download header files
    download("md5.check")
    download("md5.txt")
    download("readme.txt")
    # use file list from md5.txt to download full dataset
    with open('md5.txt', 'r') as f:
        for line in f:
            line = line.split("  ")
            file = line[1].strip()
            download(file)
