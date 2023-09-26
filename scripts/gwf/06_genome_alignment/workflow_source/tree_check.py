#!/bin/env python

import sys

newick_tree = sys.argv[1]



class TreeNode(object):
    def __init__(self, length = None, name = None, children = None):
        self.length = length
        self.name = name
        self.children = children

