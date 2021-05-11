/**
 * @file ProblemaC.cpp
 *
 * @author Ana Raquel Ferreira
 * Contact: work.raquelferreira@gmail.com
 *
 */


#include <iostream>
#include<vector>
#include<iterator>
#include<algorithm>
#include <stdio.h>
#include <queue>

using namespace std;
const int maxEquip = 1001;
#define INF 0x3f3f3f3f  // infinity

// fully connected network: ALL servers are connected: {1, 4, 7} -> {1->4, 1->7, 4->7}

// to store edge information
class Edge {
public:
    int cable;
    pair<int, int> nodes; // first->source, second->destination

    Edge(int source, int dest, int cable) {
        this->cable = cable;
        this->nodes = make_pair(source, dest);
    }
    
    // default constructor
    Edge() {
        this->cable = -1;
        this->nodes = make_pair(-1, -1);
    }

    // used to compare Edges
    bool operator==(const Edge& edge) {
        int cable = edge.cable;
        int source = edge.nodes.first;
        int dest = edge.nodes.second;
        return cable == this->cable && source == this->nodes.first && dest == this->nodes.second;
    }
};

// information for shortest path
class Spath {
public:
    vector<int> parent;
    vector<int> distance;

    Spath(vector<int> parent, vector<int> distance) {
        this->parent = parent;
        this->distance = distance;
    }

    Spath() {
        this->parent = vector<int>();
        this->distance = vector<int>();
    }
};

// A structure to represent a subset for union-find  
class subset {
public:
    int parent;
    int rank;
};

// global variables for simplicity
int servers, cableNetwork, cableTree;
vector<int> artPoints = vector<int>();

vector<Edge> edgeList = vector<Edge>();
vector<Edge> edgeAux = vector<Edge>();
vector<Edge> treeEdges = vector<Edge>();


// -------------------- aux functions --------------------

// clear adj list information
void clearAuxList(vector< pair<int, int> > auxList[]) {
    for (int u = 1; u < 1001; u++) {
        if (!auxList[u].empty()) {
            auxList[u].clear();
        }
    }
}

// aux function to return the minimum of two integers
int min(int a, int b) {
    if (a < b) {
        return a;
    }
    else {
        return b;
    }
}


// -------------------- edge functions --------------------

// add vertex to list of aux edges
void add2EdgeAux(int vertex) {
    int source, dest;
    for (int i = 0; i < (int)edgeList.size(); ++i) {
        source = edgeList[i].nodes.first;
        dest = edgeList[i].nodes.second;

        if (vertex == source || vertex == dest) {
            //printf("vertex %d found and added to edge list\n", vertex);
            // if edge has not been added, but should
            if (find(edgeAux.begin(), edgeAux.end(), edgeList[i]) == edgeAux.end()){ 
                edgeAux.push_back(edgeList[i]);
            }
        }
    }
}

// add created edge to list of tree edges
void add2TreeEdges(Edge edge) {

    //printf("edge here: source %d, dest %d, cable %d\n", edge.nodes.first, edge.nodes.second, edge.cable);
    // if edge has not been added, but should
    if (find(treeEdges.begin(), treeEdges.end(), edge) == treeEdges.end()) {
        //printf("adding?\n");
        treeEdges.push_back(edge);
    }
}

// add nodeB into the list nodeA, and nodeA into list nodeB
void add_edge(vector< pair<int, int> > list[], int nodeA, int nodeB, int cable) {
    list[nodeA].push_back(make_pair(nodeB, cable));
    list[nodeB].push_back(make_pair(nodeA, cable));
}


// -------------------------- Kruskal -------------------------- 

// A utility function to find set of an element i (uses path compression technique)  
int find(subset subsets[], int i)
{
    // find root and make root as parent of i  
    // (path compression)  
    if (subsets[i].parent != i) {
        subsets[i].parent = find(subsets, subsets[i].parent);
    }
    return subsets[i].parent;
}

// A function that does union of two sets of x and y (uses union by rank)   (union+find in ppts)
void Union(subset subsets[], int x, int y)
{
    int xroot = find(subsets, x);
    int yroot = find(subsets, y);

    // Attach smaller rank tree under root of high  
    // rank tree (Union by Rank)  
    if (subsets[xroot].rank < subsets[yroot].rank) {
        subsets[xroot].parent = yroot;
    }
    else if (subsets[xroot].rank > subsets[yroot].rank) {
        subsets[yroot].parent = xroot;
    }

    // If ranks are same, then make one as root and  
    // increment its rank by one  
    else
    {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Compare two edges according to their weights.  
// Used in sort() for sorting an array of edges  
bool comparison(Edge a, Edge b)
{
    return a.cable < b.cable;
}

// Kruskal's algorithm to find the Minimum Spanning Tree
vector<Edge> KruskalMST(vector<Edge> edgeList, int nVertex) {

    // Sort all the edges in non - decreasing order of distance
    sort(edgeList.begin(), edgeList.end(), comparison);

    int E = edgeList.size();
    
    // initialize mst
    vector<Edge> mst(nVertex);

    // Allocate memory for creating V subsets  
    subset* subsets = new subset[(nVertex * sizeof(subset))];
    
    int e = 0; // An index variable, used for mst
    int i = 0; // An index variable, used for sorted edges


    // Create V subsets with single elements
    for (int v = 0; v < nVertex; ++v)
    {
        subsets[v].parent = v;
        subsets[v].rank = 0;
    }
    
    // Number of edges to be taken is equal to V-1
    while (e < nVertex - 1 && i < E)
    {
        // Step 2: Pick the smallest edge. And increment the index for next iteration
        Edge next_edge = edgeList[i++];

        int x = find(subsets, next_edge.nodes.first);       // source
        int y = find(subsets, next_edge.nodes.second);      // dest

        // If including this edge does't cause cycle, include it in result and increment the index of result for next edge
        if (x != y)
        {
            mst[e++] = next_edge;
            Union(subsets, x, y);
        }
        // Else discard the next edge
    }

    delete[] subsets;
    return mst;
}

// -------------------- Dijkstra --------------------

// save path from one node to another (branch issue solution)
vector<int> path(vector<int> parent, int dest, vector<int> caminho) {
    // Base Case : If dest is source
    if (parent[dest] == -1) { 
        //printf("%d\t", dest);
        caminho.push_back(dest);
        return caminho; 
    }

    caminho = path(parent, parent[dest], caminho);
    caminho.push_back(dest);
    //printf("%d\t", dest);
    return caminho;
}

// Dijkstra
Spath shortestPath(vector<pair<int, int> > adj[], int nVertex, int src) {

    // create a priority queue
    priority_queue< pair<int, int>, vector <pair<int, int>>, greater<pair<int, int>> > pq;

    // vector to store distances [size of number of vertices]
    vector<int> dist(nVertex, INF);

    // insert source in pq and initialize distance to itself as zero
    pq.push(make_pair(0, src));
    dist[src] = 0;

    // vector to check if node has been analized 
    vector<bool> f(nVertex, false);

    // vector to keep the "parent" of the node (previous node visited in shortest path)
    vector<int> parent(nVertex, -1);

    // Looping till priority queue becomes empty (or all distances are not finalized)
    while (!pq.empty()) {

        // The first vertex in pair is the minimum distance vertex, extract it from priority queue. 
        int u = pq.top().second;
        pq.pop();
        f[u] = true;

        // 'i' is used to get all adjacent vertices of a vertex
        vector< pair<int, int> >::iterator i;

        for (i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            // Get vertex label and weight of current adjacent of u.
            int v = (*i).first;
            int weight = (*i).second;

            // If there is shorted path to v through u.
            if (f[v] == false && dist[v] > dist[u] + weight)
            {
                // Updating distance of v
                dist[v] = dist[u] + weight;
                parent[v] = u;
                pq.push(make_pair(dist[v], v));
            }
        }

    }

    // Print shortest distances stored in dist[]
    /*printf("Vertex Distance from Source\n");
    for (int i = 0; i < nVertex; ++i)
        printf("%d \t\t %d\n", i, dist[i]);

    printf("shortest path to 9:\n");
    path(parent, 9);*/
    Spath info = Spath(parent, dist);
    return info;
}

// ----------------------- cable cost ---------------

void costNetworkEdges(vector<Edge> edges, vector<int> serverI, int nVertex) {
    int i, j, k, e;
    int source, dest, cable;

    // turn edge list to adjacency list
    vector< pair<int, int> > adjList[maxEquip];
    for (i = 0; i < (int)edges.size(); i++) {
        source = edges[i].nodes.first;
        dest = edges[i].nodes.second;
        cable = edges[i].cable;
        if (source != -1 && dest != -1 && cable != -1) {
            add_edge(adjList, source, dest, cable);
        }
    }

    // information about shortest path
    Spath info;

    // path from one server to another
    vector<int> caminho;

    // find shortest paths to connect servers
    for (i = 0; i < (int)serverI.size() - 1; i++) {

        // find shortest path from this server to every other node
        info = shortestPath(adjList, nVertex, serverI[i]);

        // from source to all other servers
        for (j = i + 1; j < (int)serverI.size(); j++) {
            
            // add connection to the list
            add2TreeEdges(Edge(serverI[i], serverI[j], info.distance[serverI[j]]));

            // calculate path from one server to the other
            caminho = path(info.parent, serverI[j], caminho);

            // cycle through path
            for (k = 0; k < (int)caminho.size() - 1; k++) {

                source = caminho[k];
                dest = caminho[k + 1];

                // find edges in path
                for (e = 0; e < (int)edges.size(); e++) {
                    
                    // found correct edge, add to the cable
                    if ((edges[e].nodes.first == source && edges[e].nodes.second == dest) || (edges[e].nodes.first == dest && edges[e].nodes.second == source)) {

                        // if edge is available
                        if (edges[e].nodes.first != -1 && edges[e].nodes.second != -1) {
                            cableNetwork += edges[e].cable;
                        }
                    }
                }
            }
            caminho.clear();
        }
    }
}

void costTreeEdges(vector<Edge> edges) {
    int i;
    for (i = 0; i < (int)edges.size(); i++) {
        if (edges[i].cable != -1) {
            cableTree += edges[i].cable;
        }
    }
}

// -------------------- DFS - AP -------------------- 

// add to vector of AP
void foundAP(int vertex) {
    // if vertex is not yet in the vector, add it
    if (!(find(artPoints.begin(), artPoints.end(), vertex) != artPoints.end())) {
        artPoints.push_back(vertex);
        servers++;
    }
}

void DFSAP(int vertex, bool visited[], vector< pair<int, int> > adjList[], int dfs[], int low[], int parent[], vector< pair<int, int> > auxList[]) {
    // discovery time
    static int time = 0;
    
    // mark the current node as visited
    visited[vertex] = true;

    // add the node list to the aux vector -> connected components
    auxList[vertex] = adjList[vertex];

    // add edges to aux edge list
    add2EdgeAux(vertex);

    //printf("vertex %d in dfs\n", vertex);

    // initialize discovery time and low value
    dfs[vertex] = low[vertex] = ++time;

    // count the children in DFS tree
    int children = 0;

    // go through all vertices adjacent to this one
    vector< pair<int, int> >::iterator i;
    for (i = adjList[vertex].begin(); i != adjList[vertex].end(); ++i) {

        int v = i->first; // v is current adjacent of vertex
        
        // if v is not visited yet, make it a child of vertex in DFS tree and recur
        if (!visited[v]) {
            children++;
            parent[v] = vertex;
            DFSAP(v, visited, adjList, dfs, low, parent, auxList);

            // check if subtree rooted with v has a connection to one of ancestors of vertex
            low[vertex] = min(low[vertex], low[v]);

            // vertex is AP in cases:

            // (1) -> vertex is root of DFS tree and has two or more children
            if (parent[vertex] == -1 && children > 1) {
                foundAP(vertex);
            }

            // (2) -> if vertex is not root and low value of one of its children is more than dfs value (discovery value) of vertex
            if (parent[vertex] != -1 && low[v] >= dfs[vertex]) {
                foundAP(vertex);
            }
        }
        else if (v != parent[vertex]) {
            low[vertex] = min(low[vertex], dfs[v]);
        }
    }

}

void findAP(vector< pair<int, int> > adjList[], int nVertex) {
    bool *visited = new bool[nVertex];
    int *dfs = new int[nVertex];
    int *low = new int[nVertex];
    int* parent = new int[nVertex];

    // edge list for mst
    vector<Edge> kruskalEdges;

    for (int i = 0; i < nVertex; i++) {
        visited[i] = false;
        parent[i] = -1;
    }

    // store the adjency list for each connected graph
    vector< pair<int, int> > auxList[maxEquip];

    for (int v = 1; v < nVertex; v++) {
        if (!visited[v]) {
            DFSAP(v, visited, adjList, dfs, low, parent, auxList);

            // only run analysis if there's more than one server
            if (artPoints.size() > 1) {
                // cost for fully connected network [dijkstra]
                costNetworkEdges(edgeAux, artPoints, nVertex);

                // minimum spanning tree [kruskal]                
                kruskalEdges = KruskalMST(treeEdges, nVertex);
                costTreeEdges(kruskalEdges);
            }

            // clear info about a connected component
            artPoints.clear();
            clearAuxList(auxList);
            edgeAux.clear();
            treeEdges.clear();
            kruskalEdges.clear();      
        }
    }

    delete[] visited;
    delete[] dfs;
    delete[] low;
}


// -------------------- Input --------------------
int main() {

    int input;
    while (cin >> input) {
        // read the number of network equipment
        int nEquipment = input + 1; // equipment starts at 1, not 0
        int nodeA, nodeB, cable;

        //create an array of lists whose size is maxEquipment (1000 + 1) [adjacency list]
        vector< pair<int, int> > adjList[maxEquip];

        // test variables
        servers = 0;
        cableNetwork = 0;
        cableTree = 0;

        edgeList.clear();
        
        // read the network information
        while(true) {
            cin >> nodeA;
            
            // run the analysis
            if (nodeA == 0) {
                //printGraph(adjList, nEquipment);    // debug

                findAP(adjList, nEquipment);

                if (servers == 0) {
                    printf("no server\n");
                }
                else {
                    printf("%d %d %d\n", servers, cableNetwork, cableTree);
                }
                break;
            } 
            
            // not the end of input, so add the edge
            else{
                cin >> nodeB;
                cin >> cable;
                
                // add edge
                add_edge(adjList, nodeA, nodeB, cable);

                // add edge to edge list
                edgeList.push_back(Edge(nodeA, nodeB, cable));
               
            }
        }
    }
    return 0;
}


