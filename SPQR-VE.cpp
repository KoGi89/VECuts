#include <stdio.h>
#include <stdlib.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/decomposition/StaticSPQRTree.h>
#include <chrono>

using namespace ogdf;

int main(int n_args, char** args)
{
    if (n_args!= 3) 
    {
       printf("\n usage: %s <input file> <output file>\n\n", args[0]);
       exit(-1);
    }

   /* read graph from file */
   FILE* fp = fopen(args[1],"r");
   int n,m;
   fscanf(fp,"%d %d",&n,&m);
   Graph G;
   node* vertex = new node[n];
   for(int i=0;i<n;i++){vertex[i]=G.newNode();} //create vertices
   for(int i=0;i<m;i++)
   {
      int x,y;
      fscanf(fp,"%d %d",&x,&y);
      G.newEdge(vertex[x-1],vertex[y-1]); //create edges
   }
   fclose(fp);


     using namespace std::chrono;
     high_resolution_clock::time_point t1 = high_resolution_clock::now();
     

   /* initialize all count(v) to 0 */
   int* count = new int[n];
   for(int i=0;i<n;i++){count[i]=0;}

   /* construct the SPQR tree, and collect the S-nodes */
   StaticSPQRTree T = StaticSPQRTree(G); 
   List<node> SNodes = T.nodesOfType(SPQRTree::NodeType::SNode);
  
   List<node>::iterator i;

   for(i=SNodes.begin();i!=SNodes.end();i++)
   {
      Graph cycle = T.skeleton(*i).getGraph(); //current cycle; skeleton of an S-node

      
      int length = cycle.numberOfNodes();

      /* 
         find the number of real edges;
         if v is on the cycle, then count(v) += n_realEdges - #{real edges that v is incident to} 
      */
      edge e = cycle.firstEdge();
      int n_realEdges=0;
      for(int j=0;j<length;j++)
      {
         char isVirtual = T.skeleton(*i).isVirtual(e);
         n_realEdges += !isVirtual;
         node tG = T.skeleton(*i).original(e->source());
         node sG = T.skeleton(*i).original(e->target());
         if(!isVirtual)
         {
            count[tG->index()]--;
            count[sG->index()]--;
         }
         e=e->succ();
      }      

      /* do count(v) += n_realEdges, for every v on the cycle */
      node v = cycle.firstNode();
      for(int j=0;j<length;j++)
      {
         node vG = T.skeleton(*i).original(v);
         count[vG->index()]+=n_realEdges;
         v=v->succ(); 
      }       
   }


     high_resolution_clock::time_point t2 = high_resolution_clock::now();
     duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   fp = fopen(args[2],"w");
   for(int i=0;i<n;i++){fprintf(fp,"%d\n",count[i]);}
   fprintf(fp,"%lf\n",(double)time_span.count());
   fclose(fp);
 
   printf("%f\n",time_span.count());
   return 0;	
}


















