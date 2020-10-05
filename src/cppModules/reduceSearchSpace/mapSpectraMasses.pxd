from libcpp.string cimport string 
from libcpp.vector cimport vector
from cpython cimport array
from libcpp.unordered_map cimport unordered_map

cdef extern from "mapSpectraMasses.hpp":
    cdef cppclass mappings:
        mappings() except + # constructor with no errors
        unordered_map[string, vector[string]] matchedBMasses
        unordered_map[string, vector[string]] matchedYMasses
        unordered_map[string, vector[string]] kmerToProts

    cdef cppclass protein:
        protein(string, string) except +
        string name
        string sequence

    cdef cppclass boundary:
        boundary(float, float) except +
        float lowerBound
        float upperBound

    mappings * mapBoundaries(vector[boundary], vector[protein], int)