from __future__ import annotations
from typing import Any, Iterable

class TreeNode:
    '''
    Node that contains Protein Sequence information
    '''
    
    def __init__(self, key, value):
        self.keys = [key]
        self.value = value
        self.children = []
        
    def add_child(self, key: str, value: Any) -> TreeNode:
        '''
        Add a child to the node being called on. Child will have the 
        key, value that are passed in. 

        Inputs:
            key:    (str) the key that will be returned when searching for sequences
            value:  (Any) the value that is being added to the child

        Outputs:
            (TreeNode) the new child
        '''
        new_child = TreeNode(key, value)
        self.children.append(new_child)
        self.add_key(key) 
        return new_child
        
    def has_child(self, value: Any) -> bool:
        '''
        Search the current node for a child that has exactly the input value

        Inputs:
            value:  (Any) the value to search for
        Outputs:
            (bool) True if ANY child node has the search value, False otherwise
        '''
        return any([child.value == value for child in self.children])
    
    def get_child(self, value: Any) -> TreeNode:
        '''
        Search the current node for a child that has exactly the input value and return
        the node that does contain that value

        Inputs:
            value:  (Any) the value to search for
        Outputs:
            (TreeNode or None) the child node with the input value, None if the value cannot be found
        '''
        if not self.has_child(value):
            return None
        
        for c in self.children:
            if c.value == value:
                return c
    
    def add_key(self, key: str) -> None:
        '''
        Add a key to the current node

        Inputs:
            key:    (str) A new key to add to the current node
        Outputs:
            None
        '''
        key not in self.keys and self.keys.append(key)
        
    def show(self, spaces=0) -> None:
        '''
        Recursively show this node and all its children nodes

        Inputs:
        kwargs:
            spaces:     (int) the number of spaces to prefix before showing the value
        Outputs:
            None
        '''
        pspaces = ''.join([' ' for _ in range(spaces)])
        print(f'{pspaces}|---> keys: {self.keys}, value: {self.value}')
        [c.show(spaces+1) for c in self.children]
        
class Tree:
    '''
    A datastructure intended to be used like a prefix tree
    '''
    
    def __init__(self, da_tolerance=.01):
        self.root = TreeNode(None, [])
        self.da_tolerance = da_tolerance
        
    def insert(self, key: str, sequence: Iterable) -> None:
        '''
        Add a seqeunce identified by a key to the tree

        Inputs:
            key:        (str) the key that identifies this sequence
            sequence:   (Iterable) any iterable

        Outputs: 
            None
        '''
        current_node = self.root
        
        for value in sequence:
            
            # if i get a child node back, make it current node and 
            # add a key 
            if current_node.get_child(value) is not None:
                current_node = current_node.get_child(value)
                current_node.add_key(key)
                
            # add a new child to current node
            # add child adds the key to the current node and creates the new child
            # that is returned
            else:
                current_node = current_node.add_child(key, value)
            
    def search(self, sequence: Iterable) -> list:
        '''
        Search the tree for a sequence

        Inputs:
            sequence:   (Iterable) the iterable to search for
        Outputs:
            (list) all keys associated with the sequence. If the sequence is not found, 
                    an empty list is returned
        '''
        current_node = self.root
        
        for value in sequence:
            if not current_node.has_child(value):
                return []
            
            current_node = current_node.get_child(value)
            
        return current_node.keys

    def get_children_of(self, sequence: Iterable) -> list:
        '''
        Search the tree for a sequence and return a list of values
        that, when added to the input sequnence, are found in the dictionary 
        stored in the tree

        Inputs:
            sequence:   (Iterable) the sequence to look for following values
        Outputs:
            (list) values that follow the sequence
        '''
        current_node = self.root 

        for value in sequence:
            if not current_node.has_child(value):
                return []
            
            current_node = current_node.get_child(value)

        return [c.value for c in current_node.children]
    
    
    def show(self):
        '''
        Print the tree to console
        '''
        print('root')
        [c.show() for c in self.root.children]
        
                
        