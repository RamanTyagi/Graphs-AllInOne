//
//  main.cpp
//  Graphs
//
//  Created by RAMAN TYAGI on 28/08/21.
//  Copyright © 2021 RAMAN TYAGI. All rights reserved.
//
/*
 #################THEORY : GRAPHS##########################
 ->What is Vector of Pairs?
 ->A pair is a container which stores two values mapped to each other, and a vector containing multiple number of such pairs is called a vector of pairs.
      --->  vector< pair <int,int> > vect;
            * Phla Element : vect[i].first
            * Dusra Element: vect[i].second
 
 **********************GRAPHS******************************
 USE:
 ->Cities and Roads wale questions!
 ->Lan wire bichane wale questions ya koi network based questions!
 ->Source se destination jana hai (0 -> 6 )
 *Smallest Rasta btao edges ki terms mein ? USE BFS
 
 #SHORTEST PATH IN DAG APPROACH or SIMPLY DO BFS:
 1.FIND TOPO SORT USING DFS
 2.DISTANCE ARRAY BNAO
 3.STACK SE ELEMENT POP KRTE JAO AUR DISTANCE CHECK KRTE JAO USKE CHILDREN KA USING ADJ_LIST !!!
 
 *Smallest Rasta btao weight ki terms mein ? Dijstra's
 *Sari cities mein taar bichana hai , Sari cities bhi cover hojae or taar bhi km lge?
 ---->Question of Minimum spanning trees
 ----->Prim's & kruskal's
 ->Files pdi hai kaafi sari , Ek file dusri pe depend krti hai.Kis order mein compile kregi?
 *Think of topological sort!(DFS)
 
 REPRESENTATION:
 1.Adjacency matrix
 -->2d matrix bna lo rows , columns dono pe nodes.
 -->fir matrix[u][v]=1 daldo agr u-v wali edge exist krti hai to!
 -->Ya matrix[u][v]=wt daldo agr Cost adjacency matrix banani hai!
 -->2d vale operations lgenge sare!
 -->Time complexity:O(n^2)
 -->Use when Vertices<10,000
 
 2.Adjacency list
 --> int arr[n] ka mtlb hai ki array h ek n size ki jisme elements int data type ko store krega!
 --> vector<int> Adj_list[n+1]; // Adjacency list
    ->Iska mtlb hai ki array ka hr ek element vector <int> data type ko store krega!
 *Adj_list[u].push_back(v)//Method to store in adj_list
 *Adj_list[v].push_back(u)//Method to store in adj_list
 -->Agr weight ke sath krani hai store to:
 vector<pair<int,int>> Adj_list[n+1] // Adj list
 *Adj_list[u].push_back({v,wt});
 *Adj_list[v].push_back({u,wt);
 -->Time Complexity : O(N+2E)
 
 #N->number of nodes
 #E->number of edges
 
 #QUESTION MEIN n,m given honge
 *n -> no of nodes hai!
 *m -> no of edges hai!
 ->fir next m lines mein pair de rkhe honge (u,v) ke
 
 ************************************************************BRIDGES IN A GRAPH / CUT EDGE ********************************************************************
 DEF :Bridges are those components in a graph whose removal increases the number of components in the graph.
 Approach:
 *DFS*
1. Take two arrays tins and low along with visited array
2. Condition for bridge : low[child]>tins[node] { Child baadme aa rha hai to cut kr skte hai!! } 
   -> Iska condition ka mtlb hai ki us node pe khi aur jgh se raasta nhi hai aane ka .....
   -> Is mein back edge ka concept lga hai , ki ksi aur node se bhi is particular node pe aa skte hai kya , to low ko update krdenge !!!!
Problem link : https://practice.geeksforgeeks.org/problems/bridge-edge-in-graph/1
************************************************************BRIDGES IN A GRAPH********************************************************************


 ************************************************************ARTICULATION POINT IN A GRAPH/CUT VERTEX*********************************************
 DEF : If removal of a vertex leads to increase in the number of components.
 Approach:
 *DFS*
1. Take two arrays tins and low along with visited array
2. Condition for bridge : low[child]>=tins[node]&&par!=-1 { Child baadme aa rha hai to cut kr skte hai!! } 
3. par == -1&&child>1 => Cut vertex bn skti hai ye bhi!!! ( Mtlb ek child ke dfs se dusre child visit nhi hore )
 // if chid == 1 -> Us parent node ko remove krne ke baad bhi number of components 1 rhege!!!!
Problem link : https://practice.geeksforgeeks.org/problems/articulation-point-1/0/
************************************************************ARTICULATION POINT IN A GRAPH/CUT VERTEX**********************************************


 ************************************************************Kosaraju's Algorithm for Strongly Connected Components (SCC)*************************
 SCC : Are those components in which if you start from any node you can reach every other node.
 Approach:
 *DFS*
1. Sort All nodes in order of finishing time. ( TOPO SORT )
2. Transpose the Graph. ( Edges Reversed )
// for(int i=0;i<V;i++)
        {
            visited[i]=0;
            for(auto edge:adj[i])
            {
                transpose[edge].push_back(i);
            }
        }
//	
3. DFS according to the finishing time.{comps++ for every unvisited node} (Stack collected from step 1)
Problem link : https://practice.geeksforgeeks.org/problems/strongly-connected-components-kosarajus-algo/1
************************************************************Kosaraju's Algorithm for Strongly Connected Components (SCC)**************************


************************************************************Bellman Ford Algorithm | Detect Negative Weight Cycle in Graphs*************************
-> Negative Weight Cycle wale mein Dijkstra fs Jaega -> Prority queue mein baar baar minimum distances aate hi rhege !!!!!
DEF : ->This algorithm helps in detecting the shortest distance when there are negative weights But one condition No negative cycles should be present.
      -> Neagative cycle hai to use detect bhi krwa dega!!
      -> Directed graph mein lgta hai ye!
      ->  Undirected ke case mein 2 edges(Undirected to directed converted) bnani pdengi!
 Approach:
 *DFS*
1. Relax all the edges N-1 Times
 if(dist[u]+wt<dist[v])
 dist[v] = dist[u]+wt
2. After relaxing (N-1) times , distance vector collects all the distances that are ultimately the shortest distances.
3. Now for checking negative cycle, Relax the edges one more time --> If distances decreases then there is a negative weight cycle.
(Bellman) -> https://practice.geeksforgeeks.org/problems/distance-from-the-source-bellman-ford-algorithm/0/?fbclid=IwAR2_lL0T84DnciLyzMTQuVTMBOi82nTWNLuXjUgahnrtBgkphKiYk6xcyJU
(Application) -> https://practice.geeksforgeeks.org/problems/negative-weight-cycle3504/1#
************************************************************Bellman Ford Algorithm | Detect Negative Weight Cycle in Graphs**************************
 
#################THEORY : GRAPHS##########################
 */

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
using namespace std;
/*
 APPROACH:
 1.DFS(preOrder traversal)
 -->Ek given destination pe m hu , apne neighbours se pucha kya tumhare paas destination tak jane ka raasta hai kya??
 -->Agr mre ksi bhi neighbour ke paas destination tk jaane ka raasta hai to iska mtlb mre paas bhi fir destination tk jaane ka raasta hoga!
      -> Bs fir return true krte jao!
 -->Binary tree mein jsa find vala sawal hai vsa hi ye hai.
****************IMP POINTS IN THIS QUESTION****************
 -->Agr src se neighbour ke paas gye or usse pucha ki tre paas raasta hai kya destination tk kaa.
     -> Kyuki hm bhi to iske neighbour h ye hmse hi boldega tu bta tu bhi to mera neighbour h.
     ->Iss vjh se ek infinite loop chljaega hm isse puchte rhenge or ye hmse puchta rhega.
 -->Ise problem ko solve krne ke liye hi ek visited array li h bool type ki jo ki btaegi ki kon konsi nodes pe se hm ho aaye hai , or jahan dobara nhi jana !
****************IMP POINTS IN THIS QUESTION****************
2.BFS(Level order traversal)
-->remove , marks* , work , add* strategy use krte chlo
-->If the removed element from the queue == destination then, return true.
 */
// Problem link : https://leetcode.com/problems/find-if-path-exists-in-graph/
bool hasPath(vector<pair<int,int>> Graph[],int src,int dest,bool visited[])
{
    cout<<src<<" ";
    if(src==dest)
        return true;
    visited[src]=true;
   // cout<<"\nvisited"<<visited[src]<<" ";
    for(auto vertex:Graph[src])
    {
        //cout<<"\nvertex"<<vertex.first<<" ";
        if(visited[vertex.first]==false)
        {
bool hasNbrPath = hasPath(Graph,vertex.first,dest,visited);
              if(hasNbrPath==true)
                  return true;
        }
    }
    return false;
}
/*
 APPROACH:DFS(preOrder traversal)
 ----->Ek string leli jisme starting mein src daldia
        Jo path cover krta chlega!
 -->Ek given destination pe m hu , apne neighbours se pucha kya tumhare paas destination tak jane ka raasta hai kya??
 -->Agr mre ksi bhi neighbour ke paas destination tk jaane ka raasta hai to iska mtlb mre paas bhi fir destination tk jaane ka raasta hoga!
      -> Hr ek step mein apne path mein vo node dalte gye jahan pr hm jaa rhe hai
 --->Agr src == dest hojata hai ksi bhi point pe to path ko print krdo!
****************IMP POINTS IN THIS QUESTION****************
 -->Agr src se neighbour ke paas gye or usse pucha ki tre paas raast hai kya destination tk kaa.
     -> Kyuki hm bhi to iske neighbour h ye hmse hi boldega tu bta tu bhi to mera neighbour h.
     ->Iss vjh se ek infinite loop chljaega hm isse puchte rhenge or ye hmse puchta rhega.
 -->Ise problem ko solve krne ke liye hi ek visited array li h bool type ki jo ki btaegi ki kon konsi nodes pe se hm ho aaye hai , or jahan dobara nhi jana !
 # Fir post order mein is visited wale raste ko dobara se open / false krdo!
    -->Agr ye ni kiya to kaafi nodes jo baaki possible raasto mein aa skte the vo bhi block hojanenge!
****************IMP POINTS IN THIS QUESTION****************
 */
void printAllPaths(vector<pair<int,int>> Graph[],int src,int dest,string path,bool visited[])
{
    if(src==dest)
        cout<<path<<"\n";
    visited[src]=true;
    for(auto vertex:Graph[src])
    {
        if(visited[vertex.first]==false)
        {
            printAllPaths(Graph,vertex.first,dest,path+" "+to_string(vertex.first),visited);
        }
    }
    visited[src]=false;
}
/*
 *********************VIMP*********************************
  APPROACH:-->DFS(preOrder traversal)
 ->static variable bna liye , graph ko traverse krte gye or un variables mein change krte gye!
 ------>Ek weight so far lelia jo hr edge ka weight ka bhi collect krta chlega!
  ----->Ek string leli jisme starting mein src daldia
         Jo path cover krta chlega!
  -->Ek given destination pe m hu , apne neighbours se pucha kya tumhare paas destination tak jane ka raasta hai kya??
  -->Agr mre ksi bhi neighbour ke paas destination tk jaane ka raasta hai to iska mtlb mre paas bhi fir destination tk jaane ka raasta hoga!
       -> Hr ek step mein apne path mein vo node dalte gye jahan pr hm jaa rhe hai
       -> Hr ek step mein apne weight so far mein us edge ka weight bhi dalte gye .
  --->Agr src == dest hojata hai ksi bhi point pe !
      To fir smallest path,largest path , ceil path to a criteria , floor path to that criteria , kth largest path nikal lia!
****************IMP POINTS IN THIS QUESTION****************
  -->Agr src se neighbour ke paas gye or usse pucha ki tre paas raast hai kya destination tk kaa.
      -> Kyuki hm bhi to iske neighbour h ye hmse hi boldega tu bta tu bhi to mera neighbour h.
      ->Iss vjh se ek infinite loop chljaega hm isse puchte rhenge or ye hmse puchta rhega.
  -->Ise problem ko solve krne ke liye hi ek visited array li h bool type ki jo ki btaegi ki kon konsi nodes pe se hm ho aaye hai , or jahan dobara nhi jana !
  # Fir post order mein is visited wale raste ko dobara se open / false krdo!
     -->Agr ye ni kiya to kaafi nodes jo baaki possible raasto mein aa skte the vo bhi block hojanenge!
****************IMP POINTS IN THIS QUESTION****************
 -->KTH LARGEST / SMALLEST KI BAAT HO TO PRIORITY QUEUE TO USE KRO!
 ------>Weak player & strong player strategy!(min heap se bekar player ko nikalte jao or ache player ko dalte jao)
 
 *********************VIMP*********************************
 */
static string spath;
static int spath_wt=INT_MAX;
static string lpath;
static int lpath_wt=INT_MIN;
static string cpath;
static int cpath_wt=INT_MAX; // CEIL : MINIMUM AMONGST LARGEST
static string fpath;
static int fpath_wt=INT_MIN; // FLOOR : MAXIMUM AMONGST SMALLEST
class Pair
{
public:
    int wsf;
    string psf;
};
class compareTo
{
public:
  bool operator()(Pair const &p1,Pair const &p2)
    {
        return p1.wsf>p2.wsf; // Chota wsf wala pair aaega to true dega! -->MIN HEAP BN JAEGI
    }
};
static priority_queue<Pair,vector<Pair>,compareTo> pq;
static int i=0;
void multiSolver(vector<pair<int,int>> Graph[],int src,int dest,bool visited[],int criteria,int k,string path,int wt_sofar)
{
    if(src==dest)
    {
       if(wt_sofar<spath_wt)
       {
           spath_wt=wt_sofar;
           spath=path;
       }
        if(wt_sofar>lpath_wt)
        {
            lpath_wt=wt_sofar;
            lpath=path;
        }
        if(wt_sofar<criteria&&wt_sofar>fpath_wt)
        {
            fpath_wt=wt_sofar;
            fpath=path;
        }
        if(wt_sofar>criteria&&wt_sofar<cpath_wt)
        {
            cpath_wt=wt_sofar;
            cpath=path;
        }
        if(i<k)
        {
            Pair temp;
            temp.wsf=wt_sofar;
            temp.psf=path;
            pq.push(temp);
            i++;
        }
        else if(i>=k)
        {
            if(wt_sofar>pq.top().wsf)
            {
                pq.pop();
                Pair t;
                t.wsf=wt_sofar;
                t.psf=path;
                pq.push(t);
            }
            i++;
        }
            
    }
    visited[src]=true;
    for(auto vertex:Graph[src])
    {
        if(visited[vertex.first]==false)
        {
            multiSolver(Graph, vertex.first, dest, visited, criteria, k, path+" "+to_string(vertex.first), wt_sofar+vertex.second);
        }
    }
    visited[src]=false;
}
/*
 ***********************VIMP*******************************
 ->Ek 2d array bna li sari arrays(components) ko dalne ke lie!
 ->Ek 1d array li jisme ek ek krke component dalte rhenge !
 Approach:
 -->Ek loop chalya visited pe , agr visited[i]==false hai sirf tbhi use getConnectedComponents ko pass kra gya!
     soruce(i) ko phle hi comp mein push krke bhejte rhe!
 -->If visited[i]==true hai -> Phle se ksi component ka part hai!
 -->getConnectedComponent mein hr source(i) pe recursive function lga dia , Jisse vo apne tree mein sare unvisited elements pe chla jae!
 -->Jse hi ke component mila use comps mein push krdia!
 
 ## IS GRAPH CONNECTED?(kya hm ek point pe khde hokr baaki sare points pe jaa skte hai?directly or indirectly)
 ---->Agr components hi 1 aaye iska mtlb hai saare connected hi h!(check for comps ka size == 1)
***********************VIMP*******************************
 */
void getConnectedComponents(vector<pair<int,int>> Graph[],int src,vector<int> &comp,bool visited[])
{
    visited[src]=true;
    for(auto vertex:Graph[src])
    {
        if(visited[vertex.first]==false)
        {
            comp.push_back(vertex.first);
            getConnectedComponents(Graph, vertex.first, comp, visited);
        }
    }
}
/*
 ************************IMP*******************************
 #Hamiltonian path -->src se start kre or saari vertices ko visit krde bina kisi vertex ko 2 baar visit kre!
 #Hamiltonian cycle -->Hamiltonian path jahan jakr ruka agr vhn se src tk ka direct rasta h , to vo path hamiltonian cycle h
 Approach:
 ->ek src se start kra DFS lgana!
 ->Ek count variable lelia check krne ke liye ki kya saari nodes visit hogyi!
 ->saari nodes jse hi visit hui(checked by count),path print krdo
 ->Ab check krna h kya sirf hamiltonian path h ya hamiltonian cycle bhi h
 ->Agr akhiri wale element se source tk ka direct raasta hai to hamiltonian cycle(loop chlakr dekhlo us last soruce pe ki kya original source usme dla hua hai?)
************************IMP*******************************
 */
// Problem link : https://practice.geeksforgeeks.org/problems/hamiltonian-path2522/1#
void hamiltonianPathAndCycle(vector<pair<int,int>> Graph[],int src,bool visited[],int cnt,int tot_vertices,string psf,int osrc)
{
    if(cnt==tot_vertices)
    {
        cout<<psf;
        for(auto e:Graph[osrc])
        {
            if(e.first==src)
            {
                cout<<"*\n";
                return;
            }
        }
        cout<<".\n";
        return;
    }
    visited[src]=true;
    for(auto vertex: Graph[src])
    {
        if(visited[vertex.first]==false)
        {
         hamiltonianPathAndCycle(Graph, vertex.first, visited, cnt+1,tot_vertices, psf+" "+to_string(vertex.first), osrc);
        }
    }
    visited[src]=false;
}
class Pair_BFS
{
public:
    int v;//vertex
    string psf;//path so far
    Pair_BFS(int vertex,string path)
    {
        v=vertex;
        psf=path;
    }
};
/*
 #####################BFS##################################
 -->Similar to level order!
 -->BFS radius mein grow krta hai!(Use queue)!
 -->DFS depth mein grow krta hai!(Use stack)!
 --->Continue is also a loop control statement just like the break statement. continue statement is opposite to that of break statement, instead of terminating the loop, it forces to execute the next iteration of the loop.
 #####################BFS##################################
 Approach:
 -->r m* w a*(remove -> mark(condition:if not visited) -> Perform work ->add(condition:if not visited)
 #Visited mein mark or use check dhyaan se krna hai!
 -->src yahan pr change nii hoga baar baar , Ye recursive function nii h!
 ** Mistake : visited[src]=true(wrong) -> visited[p.v]=true(right)
 */
void BFS(vector<pair<int,int>> Graph[],int src,bool visited[])
{
    queue<Pair_BFS> qt;
    Pair_BFS temp(src,to_string(src));
    qt.push(temp);
    while(qt.size()>0)
    {
        Pair_BFS p = qt.front();
        qt.pop();
        if(visited[p.v]==true)
        {
            continue;
        }
        visited[p.v]=true;
        cout<<p.v<<"@"<<p.psf<<"\n";
        for(auto edge:Graph[p.v])
        {
            if(visited[edge.first]==false)
            {
                Pair_BFS t(edge.first,p.psf+to_string(edge.first));
                qt.push(t);
            }
        }
    }
}
/*
 Approach:
 ##BFS approach h same!
 -->r m* w a*(remove -> mark(condition:if not visited) -> Perform work ->add(condition:if not visited)
 --> Agr ek element phle se visit ho chuka h ,queue se remove krte time.
 -->To return true krdo! nhi to return false
 -->Ek baar true return hote hi pta chl gya cycle hai to ek varible mein store kralia true , fir use print kra lenge!
 
 ####### CHECKING A CYCLE IN A DAG #########
 ** DFS : For detecting a cycle in a directed cyclic graph use two visited arrays, one for overall visited nodes and the other one for the current state visited nodes**
 	  -> Us wale DFS mein agr visited hogyi tb cycle ki baat krna!!!
 ** BFS : Using Kahn's alorithm(Topological sort using indegree array)
 */
bool checkCyclic(vector<pair<int,int>> Graph[],int src,bool visited[])
{
    queue<Pair_BFS> qt;
    Pair_BFS temp(src,to_string(src));
    qt.push(temp);
    while(qt.size()>0)
    {
        Pair_BFS p=qt.front();
        qt.pop();
        if(visited[p.v]==true)
        {
            return true;
        }
        visited[p.v]=true;
        for(auto edge:Graph[p.v])
        {
            if(visited[edge.first]==false)
            {
                Pair_BFS t(edge.first,p.psf
                           +to_string(edge.first));
            }
        }
    }
    return false;
}
class Pair_bi
{
public:
    int v;
    string psf;
    int l;
    Pair_bi(int vertex,string path,int level)
    {
        v=vertex;
        psf=path;
        l=level;
    }
};
/*
 *************************VIMP*****************************
 -->A graph is called bipartite if it's possible to split it's vertices  in two sets of mutually exclusive and exhaustive subsets such that all edges are across sets.
     (Within a set koi edge na mile!)
 -->Every non-cyclic graph is bipartite (2 sets mein aaram se daal paoge vertices.)
 -->If length of cycle is even then graph is bipartite!
 -->If length of cycle is odd then graph is not bipartite!
 Reason : Odd cycle mein jo extreme point hoga vo dono sets mein aana chahta hai (See by making a graph of pentagon shape)
 Approach:
 1.BFS Approach(This following approch or using color array with only 2 colors)
 ->Visited array int type ki li , Usme starting mein -1 dala
 ->Fir levels dalte rhenge jse jse visit krenge!
 ->Agr ek vertex phle se hi visit ho chuka hai or queue mein dobara aa rha hai --> Making a cycle
 ->Ab cycle odd hai ya even ise check krne ke liye levels bna liye!
 ->Agr jo queue se abhi nikala hai uska level or visited mein jo us vertex ka level same nhi hai iska mtlb odd length ki cycle bn rhi hai.
 2.DFS Approach
   1.Ya to graph ka size nikal lo edge+1 kr rhe hai na hr step mein vo length ke lie hi hai(par.len-node.len)
   2.Use color array of 2 colors
*************************VIMP******************************
 */
// Problem Link : https://leetcode.com/problems/is-graph-bipartite/
bool checkbiPartite(vector<pair<int,int>> Graph[],int src,int visited2[])
{
    queue<Pair_bi> qt;
    Pair_bi t(src,to_string(src),0);
    qt.push(t);
    while(qt.size()>0)
    {
        Pair_bi p=qt.front();
        qt.pop();
        if(visited2[p.v]>-1)
        {
            //work
            if(p.l!=visited2[p.v])
            {
                return false;
            }
        }
        visited2[p.v]=p.l;
        for(auto edge:Graph[p.v])
        {
            if(visited2[edge.first]==-1)
            {
                Pair_bi temp(edge.first,p.psf+to_string(edge.first),p.l+1);
                qt.push(temp);
            }
        }
    }
    return true;
}
class Pair_Si
{
public:
    int v;
    string psf;
    int t;
    int cnt;
    Pair_Si(int vertex,string path,int cur_time,int per_count)
    {
        v=vertex;
        psf=path;
        t=cur_time;
        cnt=per_count;
    }
};
/*
 **********************IMP*********************************
 Que:
 1.You are given a graph  ,representing people and their connectivity.
 2.You are also given a source perosn (Who got infected) and time 't'.
 3.You are required to find how many people will get infected in time 't', if the infection spreads to neighbours of infected person in 1 unit of time.
 Sol:
 ->Ek pair class bnakr jisme vertex,path so far,cur_time,cur_count daal paye
 Approach:
 ->remove , mark(*) , work , add(*)
 ->Visited array int type ki bnakr usme time dalte gye ki konsa element kis time visit hua h!
 ->Jse hi remove wale element ka time == to_time hojae
   jo person count kr rkhe h unhe return ya print krdo!
 **Agr visited pe -1 pda hai abhi to iska mtlb hai vo person abhi infect nii hua !
 **********************IMP*********************************
 */
 void countingInfectedPerson(vector<pair<int,int>> Graph[],int src,int tot_time,int visited3[])
{
    int counter=INT_MIN;
    queue<Pair_Si> qt;
    Pair_Si t(src,to_string(src),1,1);
    qt.push(t);
    while(qt.size()>0)
    {
        Pair_Si p=qt.front();
        qt.pop();
        if(p.t==tot_time)
        {
            counter=max(p.cnt,counter);
        }
        visited3[p.v]=p.t;
        for(auto edge : Graph[p.v])
        {
            if(visited3[edge.first]==-1)
            {
                Pair_Si temp(edge.first,p.psf+to_string(edge.first),p.t+1,p.cnt+1);
                qt.push(temp);
            }
        }
    }
    cout<<counter;
}
 class Pair_Dij
 {
 public:
     int v;//vertex
     string psf;//path so far
     int wsf;//Weight so far
     Pair_Dij(int vertex,string path,int weight)
     {
         v=vertex;
         psf=path;
         wsf=weight;
     }
 };
class compareTo2
{
public:
    bool operator()(Pair_Dij const& p1,Pair_Dij const& p2)
    {
        return p1.wsf>p2.wsf; // min heap : smaller element ko priority!
    }
};
/*
 Que:
 1.You are given a graph & a source vertex.The vertices represent cities and the edges represent distance in kms
 2.You are required to find the shortest path to reach each city(in terms of kms) from the soruce city along with the total distance on path from source to destination.
 ************************DIJKSTRA**************************
 ->Shortest Path in terms of edges (no. of edges) : Use BFS
 ->Shortest Path in terms of edge weight : Use Dijkstra
 ->Approach : remove -> add(*) -> work -> mark(*)
 ->Use priority queue(min heap)
    ->Sbse km weight wale ko nikal kr print krega!
************************DIJKSTRA**************************
 */
void Dijkstra(vector<pair<int,int>> Graph[],int src,bool visited[])
{
    priority_queue<Pair_Dij,vector<Pair_Dij>,compareTo2> pq;
    Pair_Dij t(src,to_string(src),0);
    pq.push(t);
    while(pq.size()>0)
    {
        Pair_Dij p = pq.top();
        pq.pop();
        if(visited[p.v]==true)
        {
            continue;
        }
        visited[p.v]=true;
        cout<<p.v<<"via"<<p.psf<<"@"<<p.wsf;
        cout<<"\n";
        for(auto edge : Graph[p.v])
        {
            if(visited[edge.first]==false)
            {
                Pair_Dij temp(edge.first,p.psf+to_string(edge.first),p.wsf+edge.second);
                pq.push(temp);
            }
        }
    }
}
 class Pair_MST
 {
 public:
     int v;//vertex
     int aq_v;//acquiring vertex
     int w;//Weight
     Pair_MST(int vertex,int aq_vertex,int weight)
     {
         v=vertex;
         aq_v=aq_vertex;
         w=weight;
     }
 };
class compareTo3
{
public:
    bool operator()(Pair_MST const& p1,Pair_MST const& p2)
    {
        return p1.w>p2.w;
    }
};
/*
 **********************PRIM'S ALGORITHM********************
 Question:
 1.You are given a graph and a source vertex.The vertices represent computers and the edges represent length of the lan wire reuqired to connect them.
 2.You are required to find the minimum length of wire required to connect all PC's server and network.Print the O/P in terms of which all PC's need to be connected, and the length of wire between them.
 sol:
 ##MST
 ->Sub graph
 ->Tree(connected Acyclic))
 ->Spanning(All vertices included)
 -> N vertices and N-1 edges
 -> Every node reachable from every other node!!!
 Approach:
 *remove -> mark(*) -> work -> add(*)
 Approach:
 1.
 ->Priority queue use kri hai jisse hme sbse km weight wala pair mile har baar pop krne pr
 ->Source element ke lie acquiring vertex -1 maanli!
 ->is aq_v ke lie koi kaam nhi krna !
 2.
 -> Using Key array for weights,MST array for checking whether a particular node is in our MST,parent array for storing parents.
 -> there should be n-1 edges in MST , so outer loop will be for edges . In inner loop compare for weights and store in MST.
 -> Use priority queue for optimisation
 
 **********************KRUSKAL'S ALGORITHM**********************
 1.sort according to the weight
 2.Greedily pick up the next edge
 3. Check for their parents
 4. If same parent then , continue because cycle bnjaegi nhi to!!!!!
 5. Path compression hojaege jb find function call kra dono vertices pe.
 6. Unionn call krdo dono vertices pe!!!
 Question link : https://practice.geeksforgeeks.org/problems/minimum-spanning-tree/1
 ## CODE FOR KRUSKAL'S : 
 class Pair
{
public:
        int v1;
        int v2;
        int weight;
};
class Solution
{
    static bool comp(Pair const& a,Pair const& b)
    {
        return a.weight<b.weight;
    }
	public:
	int find(int x,vector<int>& parent)
	{
	    if(x==parent[x])
	    return x;
	    return parent[x]=find(parent[x],parent);
	}
	void unionn(int x,int y,vector<int>& parent,vector<int>& rank)
	{
	    if(rank[x]>rank[y])
	    {
	        parent[y]=x;
	    }
	    else if(rank[x]<rank[y])
	    {
	        parent[x]=y;
	    }
	    else if(rank[x]==rank[y])
	    {
	        parent[x]=y;
	        rank[y]++;
	    }
	}
	//Function to find sum of weights of edges of the Minimum Spanning Tree.
    int spanningTree(int V, vector<vector<int>> adj[])
    {
        // code here
        vector<Pair> Adj_list;
       for(int i=0;i<V;i++)
       {
           for(auto edge:adj[i])
           {
               int n = edge[0];
               int w = edge[1];
               Pair p;
               p.v1 = i;
               p.v2 = n;
               p.weight = w;
               Adj_list.push_back(p);
           }
       }
       sort(Adj_list.begin(),Adj_list.end(),comp);
        vector<int> parent(V,0);
        vector<int> rank(V,0);
        for(int i=0;i<V;i++)
          parent[i]=i;
          int sum = 0;
        for(int i=0;i<Adj_list.size();i++)
        {
           int a = Adj_list[i].v1;
           int b = Adj_list[i].v2;
           int c = Adj_list[i].weight;
           int u = find(a,parent);
           int v = find(b,parent);
           if(u==v)
            continue;
            unionn(u,v,parent,rank);
            sum+=c;
        }
        return sum;
    }
};
 */
// PRIM'S KA CODE!!!!
 void MST(vector<pair<int,int>> Graph[],int src,bool visited[])
{
    priority_queue<Pair_MST,vector<Pair_MST>,compareTo3> pq;
    Pair_MST t(src,-1,0);
    pq.push(t);
    while(pq.size()>0)
    {
        Pair_MST p = pq.top();
        pq.pop();
        if(visited[p.v]==true)
            continue;
        visited[p.v]=true;
        if(p.aq_v !=-1)
        cout<<p.v<<"-"<<p.aq_v<<"@"<<p.w<<"\n";
        for(auto edge:Graph[p.v])
        {
            if(visited[edge.first]==false)
            {
                Pair_MST temp(edge.first,p.v,edge.second);
                pq.push(temp);
            }
        }
    }
}
/*
*****************TOPOLOGICAL SORT : USING DFS*************
 Question:
 1.You are given a directed acyclic graph.The vertices represent tasks and edges represent dependencies between tasks.
 2.You are required to find and print the order in which tasks could be done.The task that should be done at last should be printed first and the task which should be done first should be printed last.This is called Topological sort.
 Solution:
 Topological sort : A permutation of vertices for a directed acyclic graph is called topological sort if for all directed edges uv , 'u' appears before 'v' in the graph.
(Agr u->v pe dependent hai na to topogical sort mein 'u','v' se upr aana chaie!)
 -->Order of compilation will be opposite to topological order!
 ->Post order mein element ko stack mein daldo=>Iska mtlb ye h ki current element jis jis pr bhi depend krta hoga vo phle se stack mein push kr chuke hai!
 --->preorder mein print krne mein dikkt aajaegi , jo element ksi or pe depend krta hai vo niche aa skta hai!
 --->post order mein print krne mein bhi order ki dikkat aaegi!
*****************TOPOLOGICAL SORT : USING DFS*************
 */
void  topologicalSort(vector<pair<int,int>> Graph[],int src,bool visited[],stack<int> &st)
{
    visited[src]=true;
    for(auto edge:Graph[src])
    {
        if(visited[edge.first]==false)
        {
            topologicalSort(Graph, edge.first, visited, st);
        }
    }
    st.push(src);
}
class Pair_it
{
public:
    int v;
    string psf;
    
};
/*
 Question:
 1.You are given a graph,and a src vertex.
 2.You are required to do a iterative depth first traversal and print which vertex is reached via which path starting fromt the src.
 solution:
 -->BFS wali r->m(*)->w->a(*) wali approach hi use krni hai bs queue ki jgh stack use krlo! kyuki depth mein jana hai hme.
 -->BFS siblings ko zyda preference deta hai , DFS child ko
 -->Reverse Preorder ka euler chlega is question mein!
 -->Recursion , call stack pe nye nye functions bnati hai=> Iska chota size hota hai call stack ka or jldi bhr jaati hai!(Agr liner bda graph hua to)
 -->Bhrte hi stack overflow exception dedega!
 -->Islie iterative DFS mein bnate h kyuki jo jm khud stack bnate hai vo heap mein bnti hai!
 */
void IterativeDFS(vector<pair<int,int>> Graph[],int src,bool visited[])
{
    stack<Pair_it> st;
    Pair_it t;
    t.psf=to_string(src);
    t.v=src;
    st.push(t);
   while(st.size()>0)
   {
       Pair_it p =st.top();
       st.pop();
       if(visited[p.v]==true)
           continue;
       visited[p.v]=true;
       cout<<p.v<<"@"<<p.psf<<"\n";
       for(auto edge:Graph[p.v])
       {
           if(visited[edge.first]==false)
           {
               Pair_it temp;
               temp.v=edge.first;
               temp.psf=p.psf+to_string(edge.first);
               st.push(temp);
           }

       }
   }
}
int main() {
    // insert code here...
    int n,m;
    
    cout<<"ENTER NUMBER OF NODES:";
    cin>>n;
    cout<<"ENTER NUMBER OF EDGES:";
    cin>>m;
  //  int Adj_matrix[n+1][n+1]; // Adj matrix
  //  vector<int> Adj_list[n+1]; // Adj list without weight
    vector<pair<int,int>> Adj_list[n+1];
    for(int i=0;i<m;i++)
    {
        int u,v,w;
        cout<<"\nENTER BOTH VERTICES OF THE EDGE:";
        cin>>u;
        cin>>v;
        cout<<"\nENTER WEIGHT OF THIS EDGE:";
        cin>>w;
        //FOr undirected graph
        Adj_list[u].push_back({v,w});
        Adj_list[v].push_back({u,w});
    }
    bool visited[n];
    for(int i=0;i<n;i++)
    {
        visited[i]=false;
    }
    /*cout<<"\nIS THERE ANY PATH BETWEEN THE GIVEN SOURCE AND DESTINATION:";
    cout<<hasPath(Adj_list,0,6,visited);*/
    /*cout<<"\nPRINTING ALL THE PATHS BETWEEN THE GIVEN SOURCE AND DESTINATION:";
    printAllPaths(Adj_list,0,6,"0",visited);*/
   /* multiSolver(Adj_list,0,6,visited,42,2,"0",0);
    cout<<"\nPRINTING THE SMALLEST PATH SO FAR:";
    cout<<spath<<"@"<<spath_wt;
    cout<<"\nPRINTING THE LARGEST PATH SO FAR:";
    cout<<lpath<<"@"<<lpath_wt;
    cout<<"\nPRINTING THE PATH WHOSE WEIGHT IS CEIL TO 42";
    cout<<cpath<<"@"<<cpath_wt;
    cout<<"\nPRINTING THE PATH WHOSE WEIGHT IS FLOOR TO 42";
    cout<<fpath<<"@"<<fpath_wt;
    cout<<"\nPRINTING THE THIRD LARGEST PATH:";
    cout<<pq.top().psf<<"@"<<pq.top().wsf;*/
   /* cout<<"\nGET CONNECTED COMPONENTS OF THE GRAPH->";
    vector<vector<int>> comps;
    for(int i=0;i<n;i++)
    {
        vector<int> comp;
        if(visited[i]==false)
        {
            comp.push_back(i); getConnectedComponents(Adj_list,i,comp,visited);
            comps.push_back(comp);
        }
    }
    for(int i=0;i<comps.size();i++)
    {
        for(int j=0;j<comps[i].size();j++)
        {
            cout<<comps[i][j]<<" ";
        }
        cout<<"\n";
    }
    if(comps.size()==1)
        cout<<"\nGRAPH IS COMPLETELY CONNECTED!";*/
   // hamiltonianPathAndCycle(Adj_list,0,visited,1,7,"0",0);
    /*cout<<"\nPRINTING THE GRAPH USING BREADTH FIRST TRAVERSAL:";
    BFS(Adj_list,0,visited);*/
    /*bool isCycle=false;
    cout<<"\nCHECKING WHETHER GRAPH IS CYCLIC OR NOT!";
    for(int i=0;i<n;i++)
    {
        if(visited[i]==false)
        {
            if(checkCyclic(Adj_list,0,visited)==true)
                isCycle=true;
        }
    }
    cout<<"\nCHECKING-------->"<<isCycle;*/
   /* int visited2[n];
    for(int i=0;i<n;i++)
        visited2[i]=-1;
    cout<<"\n CHECKING WHETHER THE GRAPH IS BIPARTITE:";
    bool isbiPartite=true;
    for(int i=0;i<n;i++)
    {
        if(visited2[i]==-1)
        {
            if(checkbiPartite(Adj_list,i,visited2)==false)
            {
             isbiPartite=false;
            }
        }
    }
    cout<<"\nCHECKING-------->"<<isbiPartite;*/
   /* int visited3[n];
    for(int i=0;i<n;i++)
        visited3[i]=-1;
    cout<<"\nCOUNTING TOTAL PERSON INFECTED DUE TO THE SPREAD OF VIRUS:";
     countingInfectedPerson(Adj_list,0,3,visited3);*/
  /*  cout<<"\nSHORTEST PATH IN TERMS OF WEIGHTS IN THE GRAPH:";
    Dijkstra(Adj_list,0,visited);*/
   /* cout<<"\n CONNCETING ALL PC'S TO A NETWORK WITH MINIMUM WIRE USING PRIM'S ALGORITHMS";
    MST(Adj_list,0,visited);*/
   /* cout<<"\n ORDER OF COMPILATION:TOPOLOGICAL SORT WITH DFS:-->";
    stack<int> st;
    for(int i=0;i<n;i++)
    {
        if(visited[i]==false)
        {
            topologicalSort(Adj_list,i,visited,st);
        }
    }
    cout<<"\nPRINTING IN A TOPOLOGICAL ORDER:";
    while(st.size()>0)
    {
        cout<<st.top()<<" ";
        st.pop();
    }*/
    /*cout<<"\nPRINTING THE GRAPH USING ITERATIVE DFS:";
    IterativeDFS(Adj_list,0,visited);*/
    return 0;
}
