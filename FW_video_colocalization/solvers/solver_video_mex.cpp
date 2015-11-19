#include "mex.h"
#include "mexUtil.h"
#include <string.h>

using namespace std;


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{ 

  
  enum{ nodes_weight_i, edges_i, idframes_i };//inputs
  enum{ ass_o };//output

  myCheckArgNumber(nrhs, 3,nlhs,1);//check nb input / nb output

  /***************************************************************
  *                           INPUT
  ***************************************************************/
  
  int nNodes;// nb of nodes
  double* nodes_weight = (double*) myCheckArg(prhs, nodes_weight_i, &nNodes, 1, myDouble);

  int nEdges;// number of edges
  int* edges = (int*) myCheckArg(prhs, edges_i, &nEdges, 2, myInt);

  // id of frame containing nodes
  int* idframes = (int*) myCheckArg(prhs, idframes_i, &nNodes, 1, myInt);

  /***************************************************************
  *                           INIT SHORTEST PATH
  ***************************************************************/

  // assignement // TODO: change it to bit vector
  int*    parents = new int[nNodes];
  double* scores  = new double[nNodes];
    
  memcpy(scores, nodes_weight, sizeof(double) * nNodes); // init scores with weights

  int last_frame = idframes[nNodes-1];// last frame

  /***************************************************************
  *                           SHORTEST PATH
  ***************************************************************/

  int n, e = 0, it, p, pm, argmin = -1;
  double sm, min_score;
  bool last_frame_visited = false;
  for(n = 0; n < nNodes; n++){

    while(e < nEdges && edges[e + nEdges] < n) e++;

    bool visited = false;
    while(edges[e + nEdges] == n){
      p = edges[e];
      if( !visited || scores[p] < sm){
        pm = p;
        sm = scores[p];
        visited = true;
      }
      e++;
    }
    

    if(visited){
      parents[n] = pm;
      scores[n] += sm;
    }
    else{
      parents[n] = -1;
    }
    //printf("n = %d - score = %f - frame= %d\n", n, scores[n]);

    if( idframes[n] == last_frame  && (!last_frame_visited || scores[n] < min_score) ){
      min_score = scores[n];
      argmin = n;
      last_frame_visited = true;
    }
     
  }

  // retrieve best path:
  plhs[ass_o] = mxCreateDoubleMatrix(nNodes, 1, mxREAL);
  double* ass    =  (double*) mxGetPr(plhs[ass_o]);
  memset(ass, 0, nNodes * sizeof(double));

  n = argmin;
  it = 0;
  while(n != -1 && it < nNodes) 
  {
    //printf("n = %d - score = %f\n", n+1, scores[n]); 
    ass[n]  = 1; 
    n       = parents[n];
    it++;
  }


  delete parents;
  delete scores;
}
