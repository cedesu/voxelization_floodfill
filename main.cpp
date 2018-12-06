#include <igl/readSTL.h>
#include <igl/winding_number.h>
#include <fstream>
#include <igl/solid_angle.h>
#include <igl/mat_min.h>
#include <igl/mat_max.h>
#include <igl/read_triangle_mesh.h>
#include "math.h" 
#include <list>
#include <vector>
#include <igl/edge_lengths.h>
#include <igl/per_edge_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/signed_distance.h>

extern const int maxn=100;
extern int cnt=0;

double max(double a,double b,double c){
	double ret=a;
	if (b>ret) ret=b;
	if (c>ret) ret=c;
	return ret;
}

bool check(const Eigen::RowVector3i now, const int X, const int Y,const int Z){
	if (now(0)<0||now(0)==X||now(1)<0||now(1)==Y||now(2)<0||now(2)==Z) return false;
	else return true;
}

void scale(const Eigen::MatrixXd &V,
           const Eigen::MatrixXi &F,
		   double &resolution,
		   Eigen::RowVectorXd &min_cord,
		   int &X,
		   int &Y,
		   int &Z){
	Eigen::RowVectorXd max_cord,a;
	igl::mat_min(V,1,min_cord,a);
	igl::mat_max(V,1,max_cord,a);
	min_cord-=Eigen::RowVector3d(0.9,0.9,0.9);
	max_cord+=Eigen::RowVector3d(0.9,0.9,0.9);
	double mx=max(max_cord(0)-min_cord(0),max_cord(1)-min_cord(1),max_cord(2)-min_cord(2));
	resolution=mx/(maxn-1);
	X=floor((max_cord(0)-min_cord(0))/resolution)+1;
	Y=floor((max_cord(1)-min_cord(1))/resolution)+1;
	Z=floor((max_cord(2)-min_cord(2))/resolution)+1;
}

void fill(std::vector<double> &voxel, std::list<Eigen::RowVector3i> &list){
Eigen::MatrixXi trans(6,3);
	trans << -1,0,0,
         1,0,0,
         0,-1,0,
         0,1,0,
         0,0,-1,
         0,0,1;
    Eigen::RowVector3i now,now1,cl;
    std::list<Eigen::RowVector3i> list1;
    std::vector<Eigen::RowVector3i> closest(maxn*maxn*maxn);
    double d;
    printf("%u\n",list.size());
    while (!list.empty()){
    	now=list.front();
    	list.pop_front();
    	for (int i=0; i<6; i++){
    		now1=now+trans.row(i);
    		if (check(now1,maxn,maxn,maxn)&&voxel[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]>1000){
    			voxel[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]=0;
    			list1.push_back(now1);
    			closest[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]=now1;
			}
		}
	}
    while (!list1.empty()){
    	now=list1.front();
    	cl=closest[now(0)*maxn*maxn+now(1)*maxn+now(2)];
    	list1.pop_front();
    	for (int i=0; i<6; i++){
    		now1=now+trans.row(i);
    		if (check(now1,maxn,maxn,maxn)){
    		  d=(now1-cl).squaredNorm();
    		  d=sqrt(d);
    		  if (voxel[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]>1000)
    		    list1.push_back(now1);
			  if (voxel[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]>d){
    			voxel[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]=d;
    			closest[now1(0)*maxn*maxn+now1(1)*maxn+now1(2)]=cl;
    			//list1.push_back(std::pair(now1,now1));
			  }
			}
		}
	}
}


int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  FILE *ff=fopen("../index1.txt","r");
  char x[60];
  std::string fi="",fl="";
  for (int j=0; j<10000; j++){
  	printf("stl %u\n",j);
  	fscanf(ff,"%s",&x);
	if (x=="") break;
  	printf("%s\n",x);
  	fl=x;
  	fi="../../Thingi10K/raw_meshes/"+fl; 
  	if (j>=2800&&j<3000){
  igl::readSTL(fi,V,F,N);
  //if ((F.rows()<10000)&&!(F.rows()<5000&&V.rows()<10000)){

  if ((F.rows()<5000&&V.rows()<10000)){
  	
Eigen::MatrixXi T;
igl::AABB<Eigen::MatrixXd,3> tree;
Eigen::MatrixXd FN,VN,EN;
Eigen::MatrixXi E;
Eigen::VectorXi EMAP;
double max_distance = 1;
    Eigen::VectorXi I;
    Eigen::MatrixXd C,N;
	
	
  // Encapsulated call to point_mesh_squared_distance to determine bounds
  {
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C;
    igl::point_mesh_squared_distance(V,V,F,sqrD,I,C);
    max_distance = sqrt(sqrD.maxCoeff());
  }
  
  // Precompute signed distance AABB tree
  tree.init(V,F);
  // Precompute vertex,edge and face normals
  igl::per_face_normals(V,F,FN);
  igl::per_vertex_normals(
    V,F,igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_ANGLE,FN,VN);
  igl::per_edge_normals(
    V,F,igl::PER_EDGE_NORMALS_WEIGHTING_TYPE_UNIFORM,FN,EN,E,EMAP);
	
  
  double resolution=2;
  Eigen::RowVectorXd min_cord;
  int X,Y,Z;
  double start=time(NULL),stop;
  printf("%u %u\n",V.rows(),F.rows());
  scale(V,F,resolution,min_cord,X,Y,Z);
  stop=time(NULL);
  printf("%f\n",(double)difftime(stop,start));
  
  int t=0,tt=0;
  Eigen::MatrixXd O(X*Y*Z,3);
  Eigen::VectorXd W(X*Y*Z);
  for (int i=0; i<X; i++)
    for (int j=0; j<Y; j++)
      for (int k=0; k<Z; k++){
      	O.row(i*Y*Z+j*Z+k)=min_cord+resolution*Eigen::RowVector3d(i,j,k);
	  }
    signed_distance_pseudonormal(O,V,F,tree,FN,VN,EN,EMAP,W,I,C,N);
  std::vector<double> voxel(maxn*maxn*maxn,10000);
  std::list<Eigen::RowVector3i> list;
  int i1,j1,z1;
  for (int i=0; i<X; i++){
    for (int j=0; j<Y; j++){
	  for (int k=0; k<Z; k++)if (W(i*Y*Z+j*Z+k)<0){
	  	i1=i+(maxn-X)/2;
	  	j1=j+(maxn-Y)/2;
	  	z1=k+(maxn-Z)/2;
	    voxel[i1*maxn*maxn+j1*maxn+z1]=W(i*Y*Z+j*Z+k)/resolution;
	    list.push_back(Eigen::RowVector3i(i1,j1,z1));
	  }
	}
  }
  fill(voxel,list);
  stop=time(NULL);
  printf("step2 %f %u\n",(double)difftime(stop,start),X*Y*Z);
  std::string fo="../../Thingi10K/sdf_newnew/";
  fo=fo+fl.substr(0,fl.size()-3);
  fo=fo.substr(0,fo.size()-3)+"out";
  FILE *f=fopen(fo.c_str(),"w");
  //fprintf(f,"%d %d %d\n",X,Y,Z);
  for (int i=0; i<maxn*maxn*maxn; i++)
  	fprintf(f,"%f ",voxel[i]);
  fclose(f);
}}}

fclose(ff);
}
