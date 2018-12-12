#include <igl/opengl/glfw/Viewer.h>
#include <igl/readSTL.h>
#include <igl/winding_number.h>
#include <fstream>
#include <igl/solid_angle.h>
#include <igl/mat_min.h>
#include <igl/mat_max.h>
#include "math.h" 
#include <list>
#include <vector>

extern const int maxn=100;
extern int cnt=0;

double max(double a,double b,double c){
	double ret=a;
	if (b>ret) ret=b;
	if (c>ret) ret=c;
	return ret;
}

void scale(const Eigen::MatrixXd &V,
           const Eigen::MatrixXi &F,
		   double &resolution,
		   Eigen::RowVectorXd &min_cord,
		   int &X,
		   int &Y,
		   int &Z){
	Eigen::RowVectorXd max_cord;
	igl::mat_min(V,1,min_cord,Eigen::RowVectorXd());
	igl::mat_max(V,1,max_cord,Eigen::RowVectorXd());
	min_cord-=Eigen::RowVector3d(1,1,1);
	double mx=max(max_cord(0)-min_cord(0),max_cord(1)-min_cord(1),max_cord(2)-min_cord(2));
	if (mx/(maxn-1)>resolution) resolution=mx/(maxn-1);
	X=floor((max_cord(0)-min_cord(0))/resolution)+1;
	Y=floor((max_cord(1)-min_cord(1))/resolution)+1;
	Z=floor((max_cord(2)-min_cord(2))/resolution)+1;
}

bool check(const Eigen::RowVector3i now, const int X, const int Y,const int Z){
	if (now(0)<0||now(0)==X||now(1)<0||now(1)==Y||now(2)<0||now(2)==Z) return false;
	else return true;
}

/*void flood_fill_pre(const Eigen::MatrixXd &V,
           		const Eigen::MatrixXi &F,
		   		const double resolution,
		   		const Eigen::RowVectorXd &min_cord,
		   		const int X,Y,Z,
				int (&voxel)[X][Y][Z],
				int x,y,z){
	std::list<> lst;
    lst.push_back(Eigen::RowVectorXi(0,0,0));
    voxel=-Eigen::MatrixXi::Identity(X,Y,Z);
    voxel[0][0][0]=0;
    double w;
    while (!lst.empty()){
    	now=lst.front();
    	lst.pop_front();
    	for (int i=0; i<6; i++){
    		now+=trans.row(i);
    		if (check(now,X,Y,Z)&&voxel[now(0)][now(1)][now(2)]==-1){
    			w=igl::wingding_number(V,F,min_cord+now*resolution);
    			if (w>0.5){
					voxel[now(0)][now(1)][now(2)]=1;
					x=now(0);
					y=now(1);
					z=now(2);
					return;
				}
    			voxel[now(0)][now(1)][now(2)]=0;
    			lst.push_back(now);
			}
    		now-=trans.row(i);
		}
	}
}*/

void flood_fill_bnd(const Eigen::MatrixXd &V,
           		const Eigen::MatrixXi &F,
		   		const double resolution,
		   		const Eigen::RowVectorXd &min_cord,
		   		const int X,
				const int Y,
				const int Z,
				std::vector<int> &voxel){
	Eigen::MatrixXi trans(6,3),trans1(8,3);
	trans << -1,0,0,
         1,0,0,
         0,-1,0,
         0,1,0,
         0,0,-1,
         0,0,1;
	trans1 << 0,0,0,
          0,0,1,
          0,1,0,
          0,1,1,
          1,0,0,
          1,0,1,
          1,1,0,
          1,1,1;
    
	std::list<Eigen::RowVector3i> lst;
    double w;
    Eigen::RowVector3i now=Eigen::RowVector3i(0,0,0),now1;
    Eigen::MatrixXi O1(V.rows()*8,3);
    int o=0;
	for (int i=0; i<V.rows(); i++){
		for (int j=0; j<3; j++){
			now(j)=floor((V(i,j)-min_cord(j))/resolution);
		}
		for (int j=0; j<8; j++){
			now+=trans1.row(j);
			if (check(now,X,Y,Z)&&voxel[now(0)*Y*Z+now(1)*Z+now(2)]==-1){
				cnt++;
				//w=igl::winding_number(V,F,min_cord+resolution*Eigen::RowVector3d(now(0),now(1),now(2)));
				O1.row(o)=now;
				voxel[now(0)*Y*Z+now(1)*Z+now(2)]=2;
				o++;
			}
			now-=trans1.row(j);
		}
	}
	Eigen::MatrixXd O(o,3);
	for (int i=0; i<o; i++){
		O.row(i)=min_cord+resolution*Eigen::RowVector3d(O1(i,0),O1(i,1),O1(i,2));
	}
	Eigen::VectorXd W;
	igl::winding_number(V,F,O,W);
	int _o=0;
	for (int i=0; i<o; i++){
		now=O1.row(i);
		if (W(i)>0.5){
			voxel[now(0)*Y*Z+now(1)*Z+now(2)]=1;
		}
		else{
			voxel[now(0)*Y*Z+now(1)*Z+now(2)]=0;
		}
		lst.push_back(now);
		_o++;
	}
	bool flag;
	std::list<Eigen::RowVector3i> _lst,n1;
		printf("xex\n");
	o=_o;
    while (o>0){
    	_o=0;
    	//_lst.clear();
    	//n1.clear();
    	for (int i=0; i<o; i++){
    		now=lst.front();
    		lst.pop_front();
    		flag=false;
    		now1=now;
    		for (int j=0; j<6; j++){
    			now=now+trans.row(j);
				if (check(now,X,Y,Z)&&voxel[now(0)*Y*Z+now(1)*Z+now(2)]==-1){
					_lst.push_back(now);
					n1.push_back(now1);
					voxel[now(0)*Y*Z+now(1)*Z+now(2)]=2;
					_o++;
				}
				now=now-trans.row(j);
			}	
		}
		Eigen::MatrixXd O(_o,3);
		int i=0;
		printf("yey\n");
		for (Eigen::RowVector3i vec:_lst){
			O.row(i)=min_cord+resolution*Eigen::RowVector3d(vec(0),vec(1),vec(2));
			i++;
		}
		igl::winding_number(V,F,O,W);
		o=0;
		//lst.clear();
		Eigen::RowVector3i _now;
		for (int i=0; i<_o; i++){
			now=_lst.front();
			_lst.pop_front();
			now1=n1.front();
			n1.pop_front();
			if (W(i)>0.5){
				voxel[now(0)*Y*Z+now(1)*Z+now(2)]=1;
			}
			else{
				voxel[now(0)*Y*Z+now(1)*Z+now(2)]=0;
			}
			flag=false;
			for (int j=0; j<6; j++){
				_now=now1+trans.row(j);
				if (check(_now,X,Y,Z)&&(voxel[_now(0)*Y*Z+_now(1)*Z+_now(2)]==1||voxel[_now(0)*Y*Z+_now(1)*Z+_now(2)]==0)&&
				voxel[now1(0)*Y*Z+now1(1)*Z+now1(2)]!=voxel[_now(0)*Y*Z+_now(1)*Z+_now(2)]){
					flag=true;
				}
			}
			if (flag){
				o++;
				lst.push_back(now);
			}
		}
	}
}

void fill(const int X,
          const int Y,
		  const int Z,
		  std::vector<int> &voxel){
	double s;
	for (int i=0; i<X; i++)
		for (int j=0; j<Y; j++){
			s=0;
			for (int k=0; k<Z; k++){
				if (voxel[i*Y*Z+j*Z+k]==-1) voxel[i*Y*Z+j*Z+k]=s;
				else s=voxel[i*Y*Z+j*Z+k];
			}
		}
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,N;
  Eigen::MatrixXi F;
  igl::readSTL(
    (argc>1?argv[1]:"../../data/52134.stl"),V,F,N);
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  //viewer.data().set_mesh(V, F.row(0));
  //viewer.data().set_face_based(true);
  //viewer.launch();
  
  double resolution=2;
  Eigen::RowVectorXd min_cord;
  int X,Y,Z;
  double start=time(NULL),stop;
  printf("%u %u\n",V.rows(),F.rows());
  scale(V,F,resolution,min_cord,X,Y,Z);
  stop=time(NULL);
  printf("%f\n",(double)difftime(stop,start));
  
  std::vector<int> voxel(X*Y*Z,-1);
  flood_fill_bnd(V,F,resolution,min_cord,X,Y,Z,voxel);
  stop=time(NULL);
  printf("%f %u\n",(double)difftime(stop,start),cnt);
  fill(X,Y,Z,voxel);
  stop=time(NULL);
  printf("%f %u\n",(double)difftime(stop,start),X*Y*Z);
  int t=0,tt=0;
  for (int i=0; i<X; i++)
  	for (int j=0; j<Y; j++)
  	  for (int k=0; k<Z; k++){
  	  	if (voxel[i*Y*Z+j*Z+k]==1){
  	  	  t++;
		}
	  }
  Eigen::MatrixXd VV(t*3,3);
  Eigen::MatrixXi FF(t,3);
  for (int i=0; i<X; i++)
  	for (int j=0; j<Y; j++)
  	  for (int k=0; k<Z; k++){
  	  	if (voxel[i*Y*Z+j*Z+k]==1){
  	  	  VV.row(tt*3)=min_cord+resolution*Eigen::RowVector3d(i,j,k);
  	  	  VV.row(tt*3+1)=min_cord+resolution*Eigen::RowVector3d(i+1,j,k);
  	  	  VV.row(tt*3+2)=min_cord+resolution*Eigen::RowVector3d(i,j+1,k);
  	  	  FF(tt,0)=tt*3;
  	  	  FF(tt,1)=tt*3+1;
  	  	  FF(tt,2)=tt*3+2;
  	  	  tt++;
		}
	  }
  printf("yes%u\n",t);
  //viewer.data().add_points(VV, Eigen::RowVector3d(0, 0, 0));
  viewer.data().set_mesh(VV, FF);
  viewer.data().set_face_based(true);
  viewer.launch();
  system("pause");
}
