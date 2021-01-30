//usage: ./LinearVEogdf <input_graph> <output_file>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <ogdf/basic/Graph.h>

using namespace std::chrono;
using namespace ogdf;

int n, m;
int* edges;
Graph G;
NodeArray<List<node>> Adj;

void compute_count();
int* count;

NodeArray<int> dfs;
node* idfs;
NodeArray<List<node>> T;
NodeArray<node> p, l, low, highp, M, Mp, tempChild;
NodeArray<int>  bcount, up;
NodeArray<List<node>> invM, invMp, invHighp; 

void compute_M();
void compute_highp();
void sort_ChildrenHighp();

NodeArray<node> ufparent, representative;
NodeArray<int> ufrank;
void ufInit();
node find(node);
void unite(node,node,node);

int main(int n_args, char** args)
{

   /* read graph from file */
   FILE* fp = fopen(args[1],"r");
   fscanf(fp,"%d %d",&n,&m);
   node* vertex = new node[n];
   for(int i=0;i<n;i++){vertex[i]=G.newNode();} //create vertices
   for(int i=0;i<m;i++)
   {
      int x,y;
      fscanf(fp,"%d %d",&x,&y);
      G.newEdge(vertex[x-1],vertex[y-1]); //create edges
   }
   fclose(fp);

   high_resolution_clock::time_point t1 = high_resolution_clock::now();    

    //construct adjacency list
    Adj.init(G);
    for(edge e = G.firstEdge(); e!=nullptr; e=e->succ())
    {
       node x = e->source();
       node y = e->target();
       Adj[x].pushBack(y);
       Adj[y].pushBack(x);
    }

    compute_count();

   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

   printf("%f\n",time_span.count());

   fp = fopen(args[2],"w");
   for(int i=0;i<n;i++)
   {
      fprintf(fp,"%d\n",count[i]);
   }
   fprintf(fp,"%lf\n",(double)time_span.count());
   fclose(fp);

   return 0;
}

void compute_count()
{
   count = (int*)malloc(sizeof(int)*n);
   for(int i=0;i<n;i++){count[i]=0;}

   node r=G.firstNode();
   dfs.init(G); 
   idfs = new node[n]; 
   p.init(G);
   p[r]=nullptr;
   T.init(G);
   for(node v=G.firstNode(); v!=nullptr; v=v->succ()){dfs[v]=-1;} //initialize
   tempChild.init(G);
   l.init(G); low.init(G); bcount.init(G); up.init(G);
   for(node v=G.firstNode(); v!=nullptr; v=v->succ())
   {
      l[v]=v; low[v]=v; bcount[v]=0; up[v]=0;
   }

   node* temp_vertex = new node[n];
   ListIterator<node>* temp_out = new ListIterator<node>[n];

   int Nr=0; //current DFS label
   dfs[r]=Nr; idfs[Nr++]=r; p[r]=nullptr;
   temp_vertex[0]=r; temp_out[0]=Adj[r].begin();
   int SP=0;
   //perform DFS
   while(SP!=-1)
   {
      node v=temp_vertex[SP];
      char descend=0;
      for(ListIterator<node> i=temp_out[SP]; i!=Adj[v].end(); i++)
      {
         node u = *i;
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; p[u]=v; T[v].pushBack(u); tempChild[v]=u; 
            temp_vertex[SP+1]=u; temp_out[SP+1]=Adj[u].begin(); temp_out[SP]=i;
            descend=1; break;
         }
         if(dfs[u]<dfs[v] && u!=p[v])
         {
            if(dfs[u]<dfs[l[v]])
            {
               l[v]=u;
               if(dfs[u]<dfs[low[v]])
               {
                  low[v]=u;
               }  
            }
            bcount[v]++; up[tempChild[u]]++;
         }
         else if(v==p[u])
         {
            if(dfs[low[u]]<dfs[low[v]])
            {
               low[v]=low[u];
            }
            bcount[v]+=bcount[u];
         }
      }
      if(descend){SP++;continue;}
      bcount[v]-=up[v];
      SP--;
   }


   //gather parameters
   compute_M();
   compute_highp();

   //e is a back-edge
   for(int i=2;i<n;i++)
   {
      node v=idfs[i];
      count[p[v]->index()]+=bcount[v]==1;
   }  

   //high(u)=v
   for(node v=G.firstNode(); v!=nullptr; v=v->succ())
   {
      ListIterator<node> u = invHighp[v].begin();
      ListIterator<node> c = T[v].begin();
      while(u.valid())
      {    
         if(up[*u]!=0){u++;continue;} //check if high(u)=highp(u)
         while(!(*c==T[v].back() || dfs[*(c.succ())]>dfs[*u]))
         {
            c++;
         }
         if(low[*u]==v || dfs[*u]<=dfs[Mp[*c]] )
         {
            count[v->index()]++;
         }
         u++;
      } 
   }

   //high(u)<v
   for(node x=G.firstNode(); x!=nullptr; x=x->succ())
   {
      ListIterator<node> u = invM[x].begin();
      ListIterator<node> c = invMp[x].begin();
      while(u.valid() && c.valid())
      {
         while(c.valid() && dfs[*c]>=dfs[*u]){c++;}
         if(!c.valid()){break;}
         if(up[*u]==0 && dfs[highp[*u]]<dfs[p[*c]])
         {
            int n_edges=0;
            int h=dfs[highp[*u]];
            while(c.valid() && h<dfs[p[*c]])
            {
               while(u.valid() && dfs[*c]<dfs[*u])
               {
                  n_edges++;
                  u++;
               }
               count[p[*c]->index()]+=n_edges;
               c++;
            }
         }
         else
         {
            u++;
         }
      }
   }

   //M(u)=v
   sort_ChildrenHighp(); //sort the lists of children in decreasing order w.r.t. the highp of their elements

   for(node v=G.firstNode(); v!=nullptr; v=v->succ())
   {
      if(invM[v].size()==0){continue;}
      ListIterator<node> u = invM[v].begin();
      u++;
      ListIterator<node> c = T[v].begin();
      node min=v;
      while(u.valid() && c.valid())
      {
         min=highp[*c];
         while(u.valid() && dfs[*u]>dfs[min])
         {
            count[v->index()]++;
            u++;
         }
         min=low[*c];
         c++;
         while(c.valid() && dfs[highp[*c]]>=dfs[min])
         {
            if(dfs[low[*c]]<dfs[min]) 
            {
               min=low[*c]; 
            }
            c++;
         }
         while(u.valid() && dfs[*u]>dfs[min]){u++;}
      }
      while(u.valid())
      {
         if(dfs[*u]<=dfs[min])
         {
            count[v->index()]++;
         }
         u++;
      }
   }
   
   //M(u)>v
   for(node x=G.firstNode(); x!=nullptr; x=x->succ())
   {
      ListIterator<node> c = invMp[x].begin();
      ListIterator<node> u = invM[x].begin();
      while(c.valid() && u.valid())
      {
         while(u.valid() && dfs[*u]>=dfs[p[*c]]){u++;}
         if(!u.valid()){break;}
         if(dfs[highp[*c]]<dfs[*u])
         {
            int n_edges=0;
            ListIterator<node> first = u;
            ListIterator<node> last;
            while(u.valid() && dfs[highp[*c]]<dfs[*u])
            {
               n_edges++;
               last=u;
               u++;
            }
            count[p[*c]->index()]+=n_edges;
            c++;
            while(c.valid() && dfs[*last]<dfs[p[*c]])
            {
               while(dfs[*first]>=dfs[p[*c]])
               {
                  n_edges--;
                  first++;
               }
               count[p[*c]->index()]+=n_edges;
               c++;
            }
         }
         else
         {
            c++;
         }
      }
   }
}

//computes M, Mp, invM, and invMp
void compute_M()
{
    M.init(G);
    Mp.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){M[v]=nullptr;Mp[v]=nullptr;}
   
    NodeArray<ListIterator<node>> L;
    NodeArray<ListIterator<node>> R;
    L.init(G); R.init(G);
    NodeArray<char> setL;
    setL.init(G);
    for(node v=G.firstNode(); v!=nullptr; v=v->succ()){setL[v]=0;}

    //begin computation, in a bottom-up fashion
    //for M, ignore the root; for Mp, ignore the child of the root as well
    for(int i=n-1;i>0;i--)
    {
       node v=idfs[i]; //i = dfs[v]
      
       //fix the pointers L and R
       for(ListIterator<node> cIndx=T[v].begin(); cIndx!=T[v].end(); cIndx++)
       {
          node c = *cIndx; //c is a child of v
          if(dfs[low[c]]<i)
          {
             if(!setL[v]){L[v]=cIndx; setL[v]=1;}
             R[v]=cIndx;
          }
       }

       //find M[v]
       if(l[v]!=v){M[v]=v;}
       else if(L[v]!=R[v]){M[v]=v;}
       else
       {
          node d=M[*L[v]];
          while(1)
          {
             if(dfs[l[d]]<i){M[v]=d;break;}
             while(dfs[low[*L[d]]]>=i){L[d]++;}
             while(dfs[low[*R[d]]]>=i){R[d]--;}
             if(L[d]!=R[d]){M[v]=d;break;}
             d=M[*L[d]];
          }
       }
      
       if(i==1){break;} //v is the child of the root
      
       int dfsP = dfs[p[v]];

       //find Mp[v]
       if(dfs[l[v]]<dfsP){Mp[v]=v;continue;}
       while(dfs[low[*L[v]]]>=dfsP){L[v]++;}        
       while(dfs[low[*R[v]]]>=dfsP){R[v]--;}
       if(L[v]!=R[v]){Mp[v]=v;continue;}
       else
       {
          node d=Mp[*L[v]];
          while(1)
          {
             if(dfs[l[d]]<dfsP){Mp[v]=d;break;}
             while(dfs[low[*L[d]]]>=dfsP){L[d]++;}
             while(dfs[low[*R[d]]]>=dfsP){R[d]--;}
             if(L[d]!=R[d]){Mp[v]=d;break;}
             d=Mp[*L[d]];
          }
       }
    }

    //calculate the inverse lists, and have their elements sorted in decreasing order
    invM.init(G);
    for(int i=n-1;i>0;i--)
    {
       node v=idfs[i];
       invM[M[v]].pushBack(v);
    }   

    invMp.init(G);
    for(int i=n-1;i>1;i--)
    {
       node v=idfs[i];
       invMp[Mp[v]].pushBack(v);
    } 

    //delete[] L; delete[] R; delete[] invMnext;
}

void compute_highp()
{
   highp.init(G);
   ufInit();
   for(int i=n-1;i>=0;i--)
   {
      node v=idfs[i];
      for(ListIterator<node> j=Adj[v].begin(); j!=Adj[v].end(); j++)
      {
         node u=*j;
         if(dfs[v]<dfs[u] && v!=p[u])
         {
            node x = representative[find(u)];
            while (p[x]!=v) 
            {
               highp[x] = v;    
               node next = representative[find(p[x])];
               unite(x,p[x],next);
               x = next;
            }
         }
      }
   }

   //calculate the inverse lists invHigh, and have their elements sorted in increasing order
   invHighp.init(G);
   for(int i=2;i<n;i++)
   {
      node v=idfs[i];
      invHighp[highp[v]].pushBack(v);
   }
}

void sort_ChildrenHighp()
{
    //the inverse lists invHighp must have been computed, and sorted in increasing order
    
    for(node v=G.firstNode(); v!=nullptr; v=v->succ())
    {
       T[v].clear();
    }

    for(int i=n-1; i>-1; i--)
    {
       node x=idfs[i];
       List<node> temp;
       for(ListIterator<node> c = invHighp[x].begin(); c.valid(); c++)
       {
          temp.pushFront(*c);
       }           
       for(ListIterator<node> c = temp.begin(); c.valid(); c++)
       {
          T[p[*c]].pushBack(*c);
       }      

    }
    //delete[] TnextChild;
}

void ufInit()
{
   ufparent.init(G);
   ufrank.init(G);
   representative.init(G);
   for (node v=G.firstNode(); v!=nullptr; v=v->succ())
   {
      ufparent[v] = representative[v] = v;
      ufrank[v] = 0;
   }
}

node find(node x)
{
   node r=x;
   while(ufparent[r]!=r){r=ufparent[r];}
   while(ufparent[x]!=x){node next=ufparent[x]; ufparent[x]=r; x=next;}
   return r;
}

void unite(node x, node y, node w)
{
   node r1=find(x);
   node r2=find(y);
   if(r1==r2){representative[r1]=w;return;}
   int cmp=ufrank[r1]-ufrank[r2];
   if(cmp<0)
   {
      ufparent[r1]=r2;
      representative[r2]=w;
   }
   else if(cmp>0)
   {
      ufparent[r2]=r1;
      representative[r1]=w;
   }
   else
   {
      ufparent[r1]=r2;
      ufrank[r2]++;
      representative[r2]=w;
   }
}






