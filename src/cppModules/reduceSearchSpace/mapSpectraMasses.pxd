from libcpp.string cimport string 
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map

cdef extern from "mapSpectraMasses.hpp":
    cdef cppclass mappings:
        mappings() except + # constructor with no errors
        unordered_map[string, vector[string]] matchedBMasses
        unordered_map[string, vector[string]] matchedYMasses
        unordered_map[string, vector[string]] kmerToProts

    cdef cppclass protein:
        protein(string, string) except +
        protein() except +
        string name
        string sequence

    cdef cppclass boundary:
        boundary(float, float) except +
        boundary() except +
        float lowerBound
        float upperBound

    mappings * mapBoundaries(vector[boundary], vector[protein], int)