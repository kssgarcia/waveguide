# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 10:27:42 2025

@author: agarz
"""

import numpy as np
import os
import matplotlib.pyplot as plt


# data = np.loadtxt(full_path, skiprows=4, max_rows=68)

# plt.plot(data[:,1], data[:,2], '.-')
# plt.axis('equal')

def read_file_with_n(filename):
    result = []
    with open(filename, "r") as file:
        lines = file.readlines()
    
    i = 0
    while i < len(lines):
        items = lines[i].strip().split()
        
        if len(items) == 1:  # single item → interpret as n
            n = int(items[0])
            # print(n)
            block = []
            for j in range(i + 1, i + 1 + n):
                if j < len(lines):
                    block.append(lines[j].strip().split())
            result.append(block)
            i += n + 1  # skip past this block
        else:
            result.append(items)
            i += 1
    
    return result


# Example usage
if __name__ == "__main__":
    path = r"C:\Users\agarz\Documents\Repositorios\Ondas\BEMLIB_22.10\03_laplace\body_2d"
    
    # reading file with streamlines
    stream_file = "body_2d.str"
    stream_path = os.path.join(path, stream_file)
    data = read_file_with_n(stream_path)
    #print(data)
    
    # exclude lists with 2 or less items
    filtered_data = [item for item in data if len(item) > 2]
    
    # ploting  streamlines
    plt.figure(1)
    plt.clf()
    
    nodes_flag = True
    for item in filtered_data:
        # print(len(item))
        array = np.array(item, dtype=float)
        if nodes_flag:
            plt.plot(array[:,1], array[:,2],'.-')
            plt.axis('equal')
            nodes_flag = False
        else:
            plt.plot(array[:,1], array[:,2],'-')
            plt.axis('equal')

    # array = np.array(new_list[0], dtype=float)
    # plt.plot(array[:,1], array[:,2], '.-')
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    
    # waves_path = os.environ["WAVES_PATH"]
    # file_path = os.path.join(waves_path,"figures","stream.pdf")
    # plt.savefig(file_path, bbox_inches='tight')

    # # reading potential file
    # potential_file = "body_2d.out"
    # potential_path = os.path.join(path, potential_file)

    # data = read_file_with_n(potential_path)
    # # exclude lists with 2 or less items
    # filtered_data = [item for item in data if len(item) > 2]
