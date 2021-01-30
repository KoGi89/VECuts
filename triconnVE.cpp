#include <stdio.h>
#include <stdlib.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/Triconnectivity.h>
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

   /* get the triconnected components */
   Triconnectivity tricComp(G);

   /* prepare a stack to keep track of the vertices of the polygons (since the polygons appear as sets of edges) */
   char* found = new char[n];
   for(int i=0;i<n;i++){found[i]=0;}
   int* vertexStack = new int[n];

   /* process only the polygons; a vertex forms vertex-edge cuts with all the real edges it is not incident to */
   for(int i=0; i<tricComp.m_numComp; i++)
   {
      if(tricComp.m_component[i].m_type==Triconnectivity::CompType::polygon)
      {
         int n_realEdges=0;
         int SP=0;
         for(edge e : tricComp.m_component[i].m_edges)
         {
                int x=e->source()->index();
                int y=e->target()->index();
                if(!found[x]){found[x]=1;vertexStack[SP++]=x;}
                if(!found[y]){found[y]=1;vertexStack[SP++]=y;}
                if(tricComp.m_pGC->original(e))
                {
                   n_realEdges++; 
                   count[e->source()->index()]--; count[e->target()->index()]--;
                }
         }
         for(int t=0;t<SP;t++)
         {
            int v=vertexStack[t];
            count[v]+=n_realEdges;
            found[v]=0;
         }
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


















