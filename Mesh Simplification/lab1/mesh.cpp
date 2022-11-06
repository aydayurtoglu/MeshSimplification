#include "mesh.h"

#ifdef _WIN32
#include <Windows.h>
#include "GL\glut.h"
#define M_PI 3.141592654
#elif __APPLE__
#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#endif

#include "math.h"
#include <string>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "mesh.h"
#include <map>
#include <queue>
#include <iterator>
#include <iomanip>
#include <cmath>
#include <set>

using namespace std;

void myObjType::draw() {

    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    double longestSide = 0.0;
    for (int i = 0; i < 3; i++)
        if ((lmax[i] - lmin[i]) > longestSide)
            longestSide = (lmax[i] - lmin[i]);
    glScalef(4.0 / longestSide, 4.0 / longestSide, 4.0 / longestSide);
    glTranslated(-(lmin[0] + lmax[0]) / 2.0, -(lmin[1] + lmax[1]) / 2.0, -(lmin[2] + lmax[2]) / 2.0);
    
    if (edges)
    {
        visualizeEdges();
        glDisable(GL_LIGHTING);
        glPopMatrix();
        glEnable(GL_LIGHTING);
        return;
    }
    for (int i = 1; i <= tcount; i++)
    {
        glBegin(GL_POLYGON);
        glNormal3dv(nlist[i]);
        if (!edges)
            glColor3f(0.5, 0.5, 0.5);
        for (int j = 0; j < 3; j++){
            if (smooth)
                glNormal3dv(vnlist[tlist[i][j]]);
            
            glVertex3dv(vlist[tlist[i][j]]);
        }
        glEnd();
    
    }
    glDisable(GL_LIGHTING);
    glPopMatrix();
    glEnable(GL_LIGHTING);
}

void myObjType::writeFile(char* filename)
{
    cout << "Writing to " << filename << endl;

    ofstream myfile;
    myfile.open(filename);
    
    for (int i = 1; i <= tcount; i++){
        myfile << "f " << tlist[i][0] << ' ' << tlist[i][1] << ' ' << tlist[i][2] << "\n";
    }
    
    for (int i = 1; i <= vcount; i++){
        myfile << "v " << vlist[i][0] << ' ' << vlist[i][1] << ' ' << vlist[i][2] << "\n";
    }
    
    myfile.close();
}

void myObjType::readFile(char* filename)
{
    cout << "Opening " << filename << endl;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile.is_open()) {
        cout << "We cannot find your file " << filename << endl;
        exit(1);
    }

    string line;
    int i, j;
    bool firstVertex = 1;
    double currCood;

    tcount = 0;
    vcount = 0;
    
    while (getline(inFile, line))
    {
        if ((line[0] == 'v' || line[0] == 'f') && line[1] == ' ')
        {
            if (line[0] == 'v')
            {
                vcount++;
                i = 1;
                const char* linec = line.data();
                for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
                    while (linec[i] == ' ') i++;
                    j = i;
                    while (linec[j] != ' ') j++;
                    currCood = vlist[vcount][k] = atof(line.substr(i, j - i).c_str());
                    if (firstVertex)
                        lmin[k] = lmax[k] = currCood;
                    else {
                        if (lmin[k] > currCood)
                            lmin[k] = currCood;
                        if (lmax[k] < currCood)
                            lmax[k] = currCood;
                    }
                    i = j;
                }

                firstVertex = 0;
            }
            if (line[0] == 'f')
            {
                tcount++;
                i = 1;
                const char* linec = line.data();
                for (int k = 0; k < 3; k++) {
                    while (linec[i] == ' ') i++;
                    j = i;
                    while (linec[j] != ' ' && linec[j] != '\\') j++;
                    tlist[tcount][k] = atof(line.substr(i, j - i).c_str());
                    i = j;
                    fnlist[tcount][k] = 0;
                    while (linec[j] != ' ') j++;

                }
            }
        }
    }
    computeTriangleNormals();
    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
    fnextList();
    computeNoOfComponents();
}

void myObjType::computeTriangleNormals(){
    double u[3];
    double v[3];
    double normal[3];
    for (int i = 1; i <= tcount; i++) {
        double* p1 = vlist[tlist[i][0]];
        double* p2 = vlist[tlist[i][1]];
        double* p3 = vlist[tlist[i][2]];
        
        // calculate u and v vectors for every point
        u[0] = p2[0]-p1[0];
        u[1] = p2[1]-p1[1];
        u[2] = p2[2]-p1[2];
        
        v[0] = p3[0]-p1[0];
        v[1] = p3[1]-p1[1];
        v[2] = p3[2]-p1[2];
        
        // Take the vector cross product of u and v ( n = u x v )
        normal[0] = u[1]*v[2]-u[2]*v[1];
        normal[1] = u[2]*v[0]-u[0]*v[2];
        normal[2] = u[0]*v[1]-u[1]*v[0];
        
        nlist[i][0] = normal[0];
        nlist[i][1] = normal[1];
        nlist[i][2] = normal[2];
    }
}

void myObjType::computeStat(){
    int i;
    double minAngle = 0;
    double maxAngle = 0;
    
    for (int i = 0; i < 18; i++) {
        statMinAngle[i] = 0;
        statMaxAngle[i] = 0;
    }
    
    for (int i = 1; i <= tcount; i++) {
        double* p1 = vlist[tlist[i][0]]; // point A
        double* p2 = vlist[tlist[i][1]]; // point B
        double* p3 = vlist[tlist[i][2]]; // point C
     
        // calculating distances between points
        double sideAB = sqrt( pow(p1[0]-p2[0], 2) +  pow(p1[1]-p2[1], 2) +  pow(p1[2]-p2[2], 2) );
        double sideAC = sqrt( pow(p1[0]-p3[0], 2) +  pow(p1[1]-p3[1], 2) +  pow(p1[2]-p3[2], 2) );
        double sideBC = sqrt( pow(p2[0]-p3[0], 2) +  pow(p2[1]-p3[1], 2) +  pow(p2[2]-p3[2], 2) );
        
        // calculating angle with cosinus theorem
        double pi = 3.141592653589793238463;
        double angleA = acos( (pow(sideAB, 2) + pow(sideAC, 2) - pow(sideBC, 2)) / (2 * sideAB * sideAC) );
        angleA = angleA * (180.0 / pi); // convert radians to degrees
        
        double angleB = acos( (pow(sideAB, 2) + pow(sideBC, 2) - pow(sideAC, 2)) / (2 * sideAB * sideBC) );
        angleB = angleB * (180.0 / pi); // convert radians to degrees
        
        double angleC = 180 - angleA - angleB; // angles A+B+C = 180 for a triangle
        
        minAngle = std::min({angleA, angleB, angleC});
        maxAngle = std::max({angleA, angleB, angleC});
        
        statMinAngle[int(floor(minAngle / 10))] += 1;
        statMaxAngle[int(floor(maxAngle / 10))] += 1;
    }
    
    cout << "Statistics for Maximum Angles" << endl;
    for (i = 0; i < 18; i++)
        cout << statMaxAngle[i] << " ";
    cout << endl;
    cout << "Statistics for Minimum Angles" << endl;
    for (i = 0; i < 18; i++)
        cout << statMinAngle[i] << " ";
    cout << endl;
}

int myObjType::org(OrTri ot) {
    std::map<int, int> map = {
        {0, 0},
        {1, 1},
        {2, 2},
        {3, 1},
        {4, 2},
        {5, 0}
    };
    return tlist[idx(ot)][map[ver(ot)]];
};
int myObjType::dest(OrTri ot) { return org(sym(ot)); };

void myObjType::fnextList(){
    
    adjFaces = {}; // used in computing number of components
    
    // initialize fnlist entries with 0 first
    for ( int i = 1; i <= tcount; i++) { // i = triangle index
        for ( int j = 0; j < 3; j++){ // j = version
            fnlist[i][j] = 0;
        }
    }
    
    int v0, v1;
    std::set<int> edge;
    std::map<std::set<int>, std::set<int>> orTriMap = {};
    
    // Store the edge as a key, store its OrTri as the associated value
    for ( int i = 1; i <= tcount; i++) { // i = triangle index
        for ( int j = 0; j < 3; j++){ // j = version
            OrTri ot = makeOrTri(i,j);
            v0 = org(ot);
            v1 = dest(ot);
            
            edge = {v0, v1};
            orTriMap[edge].insert(ot);
        }
    }
    // Retrieve the adjacent triangles
    for ( int i = 1; i <= tcount; i++){
        for ( int j = 0; j < 3; j++){
            OrTri ot = makeOrTri(i,j);
            v0 = org(ot);
            v1 = dest(ot);
            
            edge = {v0, v1};
            std::set<int> adjTriangles = orTriMap[edge]; // get all the adjacent triangles

            if (adjTriangles.size() != 1) { // there is an adjacent triangle
                int tri0 = *std::next(adjTriangles.begin(), 0); // triangle at index 0
                int tri1 = *std::next(adjTriangles.begin(), 1); // triangle at index 1
                
                if ( idx(tri0) == i ) // tri0 is equal to itself
                    fnlist[i][j] = tri1;
                else
                    fnlist[i][j] = tri0;
                
                adjFaces[i].insert( idx(fnlist[i][j]) );
            }
        }
    }
}

void myObjType::visualizeEdges(){
    glDisable(GL_LIGHTING);
    for (int i = 1; i <= tcount; i++) {
        for (int j = 0; j < 3; j++){
            if (fnlist[i][j] == 0) { // we want fnlist entry to be 0
                int v0 = tlist[i][j];
                int v1 = tlist[i][(j + 1) % 3];
                glBegin(GL_LINES);
                glColor3f(1.0f, 0.0f, 0.0f); // red color
                glVertex3dv(vlist[v0]);
                glVertex3dv(vlist[v1]);
                glEnd();
            }
        }
    }
    glEnable(GL_LIGHTING);
}

void myObjType::computeNoOfComponents(){
    int numOfComponents = 0;
    std::set<int> visitedTriangles;
    
    for (int i = 1; i <= tcount; i++)
        visitedTriangles.insert(i);  // enqueue
 
    // each loop is one component
    while (!visitedTriangles.empty()) {
        std::queue<int> trianglesToTraverse;
        
        // BFS traversal
        int first = *visitedTriangles.begin();
        trianglesToTraverse.push(first);  // enqueue
        
        while (!trianglesToTraverse.empty()) {
            int id = trianglesToTraverse.front();
            trianglesToTraverse.pop();
    
            if (visitedTriangles.find(id) != visitedTriangles.end()) {
                visitedTriangles.insert(id);
                for (set<int>::iterator it = adjFaces[id].begin(); it != adjFaces[id].end(); it++) {
                    int triangle = *it;
                    trianglesToTraverse.push(triangle);
                }
            }
            visitedTriangles.erase(id);
        }
        trianglesToTraverse = {};
        numOfComponents++;
    }
    cout << "Number of components: " << numOfComponents << endl;
}

void myObjType::orientTriangles(){
    // TODO
}

void myObjType::computeVertexNormals() {
    // TODO
    smooth = !smooth;
}

void myObjType::readFileTypeSTL(char* filename)
{
    cout << "Opening " << filename << endl;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile.is_open()) {
        cout << "We cannot find your file " << filename << endl;
        exit(1);
    }

    string line;
    int i, j;
    bool firstVertex = 1;
    double currCood;

    tcount = 0;
    vcount = 0;
    
    while (getline(inFile, line))
    {
        if (line[4] == 'f' || line[12] == 'v')
        {
            if (line[12] == 'v')
            {
                vcount++;
                i = 18;
                const char* linec = line.data();
                for (int k = 0; k < 3; k++) { // k is 0,1,2 for x,y,z
                    while (linec[i] == ' ') i++;
                    j = i;
                    while (linec[j] != ' ') j++;
                    
                    std::string str(line.substr(i, j - i));
                    std::stringstream ss(str);
                    double d;
                    ss >> d;
                    currCood = vlist[vcount][k] = d;
                    
                    if (firstVertex)
                        lmin[k] = lmax[k] = currCood;
                    else {
                        if (lmin[k] > currCood)
                            lmin[k] = currCood;
                        if (lmax[k] < currCood)
                            lmax[k] = currCood;
                    }
                    i = j;
                }

                firstVertex = 0;
            }
            if (line[4] == 'f')
            {
                tcount++;
                i = 16;
                const char* linec = line.data();
                for (int k = 0; k < 3; k++) {
                    while (linec[i] == ' ') i++;
                    j = i;
                    while (linec[j] != ' ' && linec[j] != '\\') j++;
                    
                    std::string str(line.substr(i, j - i));
                    std::stringstream ss(str);
                    double d;
                    ss >> d;
                    tlist[tcount][k] = d;
                    cout << d << endl;
                    cout << tlist[tcount][k] << endl;
                    
                    i = j;
                    fnlist[tcount][k] = 0;
                    while (linec[j] != ' ') j++;

                }
            }
        }
    }
    computeTriangleNormals();
    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
    computeNoOfComponents();
}

/*
 * Big Boss Project:
 * Iterative Decimation
 */
void myObjType::simplifyMesh(double percentage) {
    // set up desired number of triangles
    int desiredTriangleCount = tcount * percentage;
    
    // remove vertices until we reach the desired triangle count
    while (tcount > desiredTriangleCount) {
        // find neighboring vertices of all the vertices
        // setup neighbors and faces arrays
        findNeighbors();
        findFaces();
        
        // compute the cost of collapse for all edges
        // setup edgeCost array
        for ( int i = 1; i <= vcount; i++ )
            findEdgeCost(i);
        
        // find the shortest edge
        int index = findEdgeWithSmallestCost();
        
        // collapse vertex to its neighbor with smallest edge
        halfCollapse(index, collapseVertex[index]);
    }
    
    fnextList();
    computeTriangleNormals();
    cout << "No. of vertices: " << vcount << endl;
    cout << "No. of triangles: " << tcount << endl;
    computeStat();
}

void myObjType::findNeighbors() {
    for (int i = 1; i <= vcount; i++)
        neighbors[i].clear();
    for (int i = 1; i <= tcount; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                if (tlist[i][j] != tlist[i][k])
                    neighbors[tlist[i][j]].insert(tlist[i][k]);
            }
        }
    }
}

void myObjType::findFaces() {
    for (int i = 1; i <= vcount; i++)
        facesOfVertices[i].clear();
    for (int i = 1; i <= tcount; i++) {
        for (int j = 0; j < 3; j++)
            facesOfVertices[tlist[i][j]].insert(i);
    }
}

void myObjType::findEdgeCost(int vertex) {
     edgeCost[vertex] = -1;
     collapseVertex[vertex] = -1;
     
     if ( neighbors[vertex].size() != 0) { // vertex has neighbors
         edgeCost[vertex] = 10000;
         
         // find the edge with the smallest cost
         for (set<int>::iterator it = neighbors[vertex].begin(); it != neighbors[vertex].end(); it++) {
             int neighbor = *it;
             double cost = findEdgeLength(vertex, neighbor);
             
             if ( cost < edgeCost[vertex] ) {
                 collapseVertex[vertex] = neighbor;
                 edgeCost[vertex] = cost;
             }
         }
     }
}

double myObjType::findEdgeLength(int vertex, int neighbor) {
    return sqrt( pow((vlist[neighbor][0] - vlist[vertex][0]), 2) +
                 pow((vlist[neighbor][1] - vlist[vertex][1]), 2) +
                 pow((vlist[neighbor][2] - vlist[vertex][2]), 2));
}

int myObjType::findEdgeWithSmallestCost() {
    int index = 1;
    for (int i = 1; i <= vcount; i++) {
        if (edgeCost[i] < edgeCost[index] )
            index = i;
    }
    return index;
}

void myObjType::halfCollapse(int vertex, int neighbor) {
    if ( neighbor != -1 ) { // vertex has neighbors
        
        std::vector<int> removedFaces;
        for (set<int>::iterator it = facesOfVertices[vertex].begin(); it!= facesOfVertices[vertex].end(); it++) {
            int face = *it;
            
            // neighbor is a vertex of the triangle
            if (tlist[face][0] == neighbor || tlist[face][1] == neighbor || tlist[face][2] == neighbor)
                removedFaces.push_back(face);
            else { // replace triangle's vertex with the neighbor of the vertex
                for (int i = 0; i < 3; i++){
                    if (tlist[face][i] == vertex){
                        tlist[face][i] = neighbor;
                        break;
                    }
                }
            }
        }
        for ( int i = 0; i < removedFaces.size(); i++) {
            int face = removedFaces.at(i);

            // remove triangle from triangle list
            // move the vertices of the triangle to the next one
            for (int j = face; j < tcount; j++) {
                for (int k = 0; k < 3; k++) {
                    tlist[j][k] = tlist[j+1][k];
                }
            }
            tcount--;
            
            // decrease the indices of the other faces to be removed
            for (int j = i + 1; j < removedFaces.size(); j++) {
                if (face < removedFaces[j])
                    removedFaces[i+1]--;
            }
        }
    }
    
    // vlist[i][j] = -1
    // remove the given vertex from vertex list
    // move the points of the deleted vertex to the next one
    for ( int i = vertex; i < vcount; i++){
        for (int j = 0; j < 3; j++)
            vlist[i][j] = vlist[i + 1][j];
    }
    vcount--;
    
    // decrease indices of vertices that are bigger than the input vertex
    for ( int i = 1; i <= tcount; i++){
        for( int j = 0; j < 3; j++){
            if ( tlist[i][j] > vertex)
                tlist[i][j]--;
        }
    }
}
