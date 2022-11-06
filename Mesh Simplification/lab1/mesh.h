#pragma once

// maximum number of vertices and triangles
#define MAXV 1000000
#define MAXT 1000000

#include <iostream>
#include <map>
#include <set>

typedef int OrTri;
typedef int tIdx;

inline OrTri makeOrTri(tIdx t, int version) { return (t << 3) | version; };
inline tIdx idx(OrTri ot) { return ot >> 3; };
inline int ver(OrTri ot) { return ot & 0b111; };

inline OrTri enext(OrTri ot) {
    std::map<int, int> enextMap = {
        {0, 1},
        {1, 2},
        {2, 0},
        {3, 5},
        {4, 3},
        {5, 4}
    };
    return makeOrTri(idx(ot), enextMap[ver(ot)]);
};
inline OrTri sym(OrTri ot) {
    std::map<int, int> symMap = {
        {0, 3},
        {1, 4},
        {2, 5},
        {3, 0},
        {4, 1},
        {5, 2}
    };
    return makeOrTri(idx(ot), symMap[ver(ot)]);
};

class myObjType {
    int vcount = 0;
    int tcount = 0;
    double vlist[MAXV][3];   // vertices list
    int tlist[MAXT][3];      // triangle list
    int fnlist[MAXT][3];     // fnext list for future (not this assignment)
    double nlist[MAXT][3];   // storing triangle normals
    double vnlist[MAXV][3];  // storing vertex normal vectors
    
    std::set<int> neighbors[MAXV];         // store neighbors of vertices
    std::set<int> facesOfVertices[MAXV];   // store the corresponding faces for each vertex
    std::map<int, std::set<int>> adjFaces; // store adjacent faces of a triangle
    double edgeCost[MAXV];                 // store cost of collapsing
    int collapseVertex[MAXV];              // store vertices to collapse to
    
    double lmax[3];          // the maximum coordinates of x,y,z
    double lmin[3];          // the minimum coordinates of x,y,z

    int statMinAngle[18]; // each bucket is  degrees has a 10 degree range from 0 to 180 degree
    int statMaxAngle[18];
    bool smooth = false;
    
public:
    
    myObjType() { vcount = 0; tcount = 0; };
    bool edges = false;
    int org(OrTri ot);
    int dest(OrTri ot);
    
    void readFile(char* filename);  // assumming file contains a manifold
    void readFileTypeSTL(char* filename); 
    void writeFile(char* filename);
    
    void draw();
    void computeStat();
    void computeTriangleNormals();
    void fnextList();
    
    void visualizeEdges();
    void computeNoOfComponents();
    void orientTriangles();
    void computeVertexNormals();
    
    void simplifyMesh(double percentage);
    void findNeighbors();
    void findFaces();
    void findEdgeCost(int vertex);
    void halfCollapse(int vertex, int neighbor);
    int findEdgeWithSmallestCost();
    double findEdgeLength(int vertex, int neighbor);
};


