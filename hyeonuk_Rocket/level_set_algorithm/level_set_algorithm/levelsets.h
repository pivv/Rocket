#include "surfaces.h"

#define MAXSIZE 1
#define MINERROR 1e-7
#define MAXSAVING 1000000


/*struct linkedlist
{
	int value;
	struct linkedlist *next;
	linkedlist(int val);
};

linkedlist::linkedlist(int val)
{
	value = val;
	next = NULL;
}

struct linkedlist *addlink(struct linkedlist *start, int newval)
{
	struct linkedlist *newlink = new struct linkedlist(newval);
	if(start==NULL) return newlink;
	struct linkedlist *point = start;
	while(point->next == NULL) point = point->next;
	point->next = newlink;
	return start;
}

struct linkedlist *deletelink(struct linkedlist *start, int delval)
{
	if(start==NULL) return NULL;
	if(start->value==delval)
	{
		delete start;
		return NULL;
	}
	struct linkedlist *point = start;
	while(point->next->value == delval)
	{
		point = point->next;
		if(point==NULL) return start;
	}
	point->next = point->next->next;
	delete point->next;
	return start;
}

struct linkedlist *deleteend(struct linkedlist *start)
{
	if(start==NULL) return NULL;
	if(start->next == NULL)
	{
		delete start;
		return NULL;
	}
	struct linkedlist *point = start;
	while(point->next->next == NULL) point = point->next;
	delete point->next;
	point->next = NULL;	
	return start;
}*/


class level_set
{
public:
	level_set(int dim, int depth, double timestep, struct surface* initsurf, struct surface* bdry);
	~level_set();

	void tree_propagate_surface_velocity(int num);
	void tree_propagate_surface();

	void reinitial_scheme();
	void iter_reinitial_scheme();

	void iter_tree_reconstruct_surface();
	
	struct surface * level_surface;
	


	//void reinitial_scheme_2d_old();
	//void iter_reinitial_scheme_2d_old();


	//double ** level_phi_2d;			// level_phi_2d[mesh_x_num+2][mesh_y_num+2]
	//double *** level_phi_3d;		// level_phi_3d[mesh_x_num+2][mesh_y_num+2][mesh_z_num+2]
	//int mesh_x_num;
	//int mesh_y_num;
	//int mesh_z_num;
	//void propagate_surface();
	//void reconstruct_surface();
	//void construct_phi();
	//void construct_phi_old();
	//void extrapolate_burning_rate();

	
	double **mat_savingtree;		// test function
	double *array_savingtree;

	int numnode;
	void savingtree();
	void savingtree_recursive(struct quadtree *current);
	void savingtree_recursive(struct octree *current);

	void signphi();

	double test_volume(struct surface * after_surface);
	void test_volume_recursive(struct quadtree *current, double *d, double *e);

    double test_l1(struct surface *after_surface);
    void test_l1_recursive(struct quadtree *current, struct surface *surf, double *d, double *e);

    double test_linf(struct surface *after_surface);
    void test_linf_recursive(struct quadtree *current, struct surface *surf, double *d, double *e);

	int link_num;
	int tree_num;
	int tree_end_num;


private:

	/* Basic Factors */

	int dimension; // 1 or 2 or 3
	double time_step;
	int accuracy_order;
	struct surface * boundary;
	int max_tree_depth;


	/* Functions for 2 dimension */

	struct pointphi_2d *tree_point_phi_2d_start;
	struct pointphi_2d *tree_point_phi_2d_tail;
	struct quadtree *tree_grid_2d;
	

	void deletelink(struct pointphi_2d *current);
	struct pointphi_2d *addlink(double x, double y, double phi, struct pointphi_2d *left, struct pointphi_2d *right, struct pointphi_2d *top, struct pointphi_2d *bottom);
	
	bool lipschitz_cond(struct quadtree *current, bool b_initial);
	void generate_tree(struct quadtree *current, bool b_initial);
	void addtree(struct quadtree *current, bool b_initial);
	void deletetree(struct quadtree *current);

	double level_computephi(double x, double y);

	void tree_reconstruct_surface_2d(struct quadtree *current, struct surface *tempsurf, int &point_num, int &edge_num);
	
	void modify_coeffs(double x1, double x2, double y1, double y2);


	/* Functions for 3 dimension */

	// 1. Tree structures

	struct pointphi_3d *tree_point_phi_3d_start;
	struct pointphi_3d *tree_point_phi_3d_tail;
	struct octree *tree_grid_3d;

	// 2. making trees

	void deletelink(struct pointphi_3d *current);
	struct pointphi_3d *addlink(double x, double y, double z, double phi, struct pointphi_3d *left, struct pointphi_3d *right, struct pointphi_3d *top, struct pointphi_3d *bottom, struct pointphi_3d *front, struct pointphi_3d *back);
	
	bool lipschitz_cond(struct octree *current, bool b_initial);
	void generate_tree(struct octree *current, bool b_initial);
	void addtree(struct octree *current, bool b_initial);
	void deletetree(struct octree *current);

	// 3. extrapolating level sets

	// 4. operators for level sets
	
	double level_computephi(double x, double y, double z);

	double tree_hg(struct pointphi_3d *point, bool pm, double si1, double si2, double si3, double sj1, double sj2, double sj3);

	// 5. propagating level sets

	// 6. reconstructing level sets

	void tree_reconstruct_surface_3d(struct octree *current, struct surface *tempsurf, int &point_num, int &edge_num);

	
	//double mesh_size;
	//double tree_hg_old(struct pointphi_2d *point, bool pm, double si1, double si2, double sj1, double sj2);

	//double ** level_burning_rate_2d;	// level_burning_rate[mesh_x_num][mesh_y_num]
	//double *** level_burning_rate_3d;		// level_burning_rate_3d[mesh_x_num][mesh_y_num][mesh_z_num]
	//int padding_num;
	//int propagate_num;
	//double g_HJ(int i, int j, int k);
	//void bounding_construct_phi(int i0);
	//bool in_bounding(int i, int j, int k);
};

void level_set::deletelink(struct pointphi_2d *current)
{
	if(current->before!=NULL) current->before->next = current->next;
	if(current->next!=NULL) current->next->before = current->before;
	else if(current->before!=NULL) tree_point_phi_2d_tail = current->before;
	else
	{
		tree_point_phi_2d_start = NULL;
		tree_point_phi_2d_tail = NULL;
	}
	if(current->left!=NULL) current->left->right = current->right;
	if(current->right!=NULL) current->right->left = current->left;
	if(current->top!=NULL) current->top->bottom = current->bottom;
	if(current->bottom!=NULL) current->bottom->top = current->top;

	delete current;

	link_num--;
}

void level_set::deletelink(struct pointphi_3d *current)
{
	if(current->before!=NULL) current->before->next = current->next;
	if(current->next!=NULL) current->next->before = current->before;
	else if(current->before!=NULL) tree_point_phi_3d_tail = current->before;
	else
	{
		tree_point_phi_3d_start = NULL;
		tree_point_phi_3d_tail = NULL;
	}
	if(current->left!=NULL) current->left->right = current->right;
	if(current->right!=NULL) current->right->left = current->left;
	if(current->top!=NULL) current->top->bottom = current->bottom;
	if(current->bottom!=NULL) current->bottom->top = current->top;
	if(current->front!=NULL) current->front->back = current->back;
	if(current->back!=NULL) current->back->front = current->front;

	delete current;

	link_num--;
}

struct pointphi_2d *level_set::addlink(double x, double y, double phi, struct pointphi_2d *left, struct pointphi_2d *right, struct pointphi_2d *top, struct pointphi_2d *bottom)
{
	struct pointphi_2d *newpoint = new struct pointphi_2d();
	if(tree_point_phi_2d_tail!=NULL) tree_point_phi_2d_tail->next = newpoint;
	newpoint->before = tree_point_phi_2d_tail;
	newpoint->left = left;
	newpoint->right = right;
	newpoint->top = top;
	newpoint->bottom = bottom;

	if(left!=NULL) left->right = newpoint;
	if(right!=NULL) right->left = newpoint;
	if(top!=NULL) top->bottom = newpoint;
	if(bottom!=NULL) bottom->top = newpoint;

	newpoint->x = x;
	newpoint->y = y;
	newpoint->phi = phi;

	if(tree_point_phi_2d_start==NULL) tree_point_phi_2d_start = newpoint;

	tree_point_phi_2d_tail = newpoint;

	link_num++;

	return newpoint;
}

struct pointphi_3d *level_set::addlink(double x, double y, double z, double phi, struct pointphi_3d *left, struct pointphi_3d *right, struct pointphi_3d *bottom, struct pointphi_3d *top, struct pointphi_3d *back, struct pointphi_3d *front)
{
	struct pointphi_3d *newpoint = new struct pointphi_3d();
	if(tree_point_phi_3d_tail!=NULL) tree_point_phi_3d_tail->next = newpoint;
	newpoint->before = tree_point_phi_3d_tail;
	newpoint->left = left;
	newpoint->right = right;
	newpoint->top = top;
	newpoint->bottom = bottom;
	newpoint->front = front;
	newpoint->back = back;

	if(left!=NULL) left->right = newpoint;
	if(right!=NULL) right->left = newpoint;
	if(top!=NULL) top->bottom = newpoint;
	if(bottom!=NULL) bottom->top = newpoint;
	if(front!=NULL) front->back = newpoint;
	if(back!=NULL) back->front = newpoint;

	newpoint->x = x;
	newpoint->y = y;
	newpoint->z = z;
	newpoint->phi = phi;

	if(tree_point_phi_3d_start==NULL) tree_point_phi_3d_start = newpoint;

	tree_point_phi_3d_tail = newpoint;

	link_num++;

	return newpoint;
}

level_set::level_set(int dim, int depth, double timestep, struct surface* initsurf, struct surface* bdry)
{
	dimension = dim;
	max_tree_depth = depth;
	time_step = timestep;
	accuracy_order = 1;
	level_surface = initsurf;
	boundary = bdry;
	
	//mesh_size = 0.1;
	//mesh_x_num = (int)floor(MAXSIZE/mesh_size);
	//mesh_y_num = (int)floor(MAXSIZE/mesh_size);
	//mesh_z_num = (int)floor(MAXSIZE/mesh_size);

	//level_phi_2d = NULL;
	//level_phi_3d = NULL;
	//level_burning_rate_2d = NULL;
	//level_burning_rate_3d = NULL;

	//padding_num = 8;
	//propagate_num = 0;

	tree_grid_2d = NULL;
	tree_point_phi_2d_start = NULL;
	tree_point_phi_2d_tail = NULL;

	tree_grid_3d = NULL;
	tree_point_phi_3d_start = NULL;
	tree_point_phi_3d_tail = NULL;


	mat_savingtree = NULL;
	array_savingtree = NULL;
	numnode = 0;			// exists for saving tree
	
	link_num = 0;
	tree_num = 1;
	tree_end_num = 1;

	if(dimension==2)
	{
		tree_grid_2d = new quadtree();
		tree_grid_2d->depth = 0;
		tree_grid_2d->length = (double)MAXSIZE;
 
		/*struct pointphi_2d *initleftbottom = addlink(0.,0.,distance_circle_point(0.,0.,0.15,0.5,0.5),NULL, NULL, NULL, NULL);
		struct pointphi_2d *initlefttop = addlink(0.,(double)MAXSIZE,distance_circle_point(0.,(double)MAXSIZE,0.15,0.5,0.5),NULL, NULL, NULL, NULL);
		struct pointphi_2d *initrightbottom = addlink((double)MAXSIZE,0.,distance_circle_point((double)MAXSIZE,0.,0.15,0.5,0.5),NULL, NULL, NULL, NULL);
		struct pointphi_2d *initrighttop = addlink((double)MAXSIZE,(double)MAXSIZE,distance_circle_point((double)MAXSIZE,(double)MAXSIZE,0.15,0.5,0.5),NULL, NULL, NULL, NULL);*/

		struct pointphi_2d *initleftbottom = addlink(0.,0.,distance_surf_point(level_surface,0.,0.,0), NULL, NULL, NULL, NULL);
		struct pointphi_2d *initlefttop = addlink(0.,(double)MAXSIZE,distance_surf_point(level_surface,0.,(double)MAXSIZE,0), NULL, NULL, NULL, NULL);
		struct pointphi_2d *initrightbottom = addlink((double)MAXSIZE,0.,distance_surf_point(level_surface,(double)MAXSIZE,0.,0), NULL, NULL, NULL, NULL);
		struct pointphi_2d *initrighttop = addlink((double)MAXSIZE,(double)MAXSIZE,distance_surf_point(level_surface,(double)MAXSIZE,(double)MAXSIZE,0), NULL, NULL, NULL, NULL);

		initleftbottom->right = initrightbottom;
		initleftbottom->top = initlefttop;

		initlefttop->right = initrighttop;
		initlefttop->bottom = initleftbottom;

		initrightbottom->left = initleftbottom;
		initrightbottom->top = initrighttop;

		initrighttop->left = initlefttop;
		initrighttop->bottom = initleftbottom;

		tree_grid_2d->phi_leftbottom = initleftbottom;
		tree_grid_2d->phi_lefttop = initlefttop;
		tree_grid_2d->phi_rightbottom = initrightbottom;
		tree_grid_2d->phi_righttop = initrighttop;

		quadtree_coeffs_saving(tree_grid_2d, initleftbottom);
		quadtree_coeffs_saving(tree_grid_2d, initlefttop);
		quadtree_coeffs_saving(tree_grid_2d, initrightbottom);
		quadtree_coeffs_saving(tree_grid_2d, initrighttop);

		generate_tree(tree_grid_2d, true);
	}

	if(dimension==3)
	{
		tree_grid_3d = new octree();
		tree_grid_3d->depth = 0;
		tree_grid_3d->length = (double)MAXSIZE;

		struct pointphi_3d *initleftbottomback = addlink(0.,0.,0.,distance_surf_point(level_surface,0.,0.,0.), NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initlefttopback = addlink(0.,(double)MAXSIZE,0.,distance_surf_point(level_surface,0.,(double)MAXSIZE,0.),NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initrightbottomback = addlink((double)MAXSIZE,0.,0.,distance_surf_point(level_surface,(double)MAXSIZE,0.,0.),NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initrighttopback = addlink((double)MAXSIZE,(double)MAXSIZE,0.,distance_surf_point(level_surface,(double)MAXSIZE,(double)MAXSIZE,0.),NULL, NULL, NULL, NULL, NULL, NULL);

		struct pointphi_3d *initleftbottomfront = addlink(0.,0.,(double)MAXSIZE,distance_surf_point(level_surface,0.,0.,(double)MAXSIZE), NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initlefttopfront = addlink(0.,(double)MAXSIZE,(double)MAXSIZE,distance_surf_point(level_surface,0.,(double)MAXSIZE,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initrightbottomfront = addlink((double)MAXSIZE,0.,(double)MAXSIZE,distance_surf_point(level_surface,(double)MAXSIZE,0.,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL);
		struct pointphi_3d *initrighttopfront = addlink((double)MAXSIZE,(double)MAXSIZE,(double)MAXSIZE,distance_surf_point(level_surface,(double)MAXSIZE,(double)MAXSIZE,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL);


		initleftbottomback->right = initrightbottomback;
		initleftbottomback->top = initlefttopback;
		initleftbottomback->front = initleftbottomfront;

		initlefttopback->right = initrighttopback;
		initlefttopback->bottom = initleftbottomback;
		initlefttopback->front = initlefttopfront;

		initrightbottomback->left = initleftbottomback;
		initrightbottomback->top = initrighttopback;
		initrightbottomback->front = initrightbottomfront;

		initrighttopback->left = initlefttopback;
		initrighttopback->bottom = initleftbottomback;
		initrighttopback->front = initrighttopfront;
		
		initleftbottomfront->right = initrightbottomfront;
		initleftbottomfront->top = initlefttopfront;
		initleftbottomfront->back = initleftbottomback;

		initlefttopfront->right = initrighttopfront;
		initlefttopfront->bottom = initleftbottomfront;
		initlefttopfront->back = initlefttopback;

		initrightbottomfront->left = initleftbottomfront;
		initrightbottomfront->top = initrighttopfront;
		initrightbottomfront->back = initrightbottomback;

		initrighttopfront->left = initlefttopfront;
		initrighttopfront->bottom = initleftbottomfront;
		initrighttopfront->back = initrighttopback;
		

		tree_grid_3d->phi_leftbottomback = initleftbottomback;
		tree_grid_3d->phi_lefttopback = initlefttopback;
		tree_grid_3d->phi_rightbottomback = initrightbottomback;
		tree_grid_3d->phi_righttopback = initrighttopback;

		tree_grid_3d->phi_leftbottomfront = initleftbottomfront;
		tree_grid_3d->phi_lefttopfront = initlefttopfront;
		tree_grid_3d->phi_rightbottomfront = initrightbottomfront;
		tree_grid_3d->phi_righttopfront = initrighttopfront;

		generate_tree(tree_grid_3d, true);
	}


	//construct_phi();

	//extrapolate_burning_rate();
}

level_set::~level_set()
{
	if(dimension == 2)
	{

		if(tree_point_phi_2d_start!=NULL)
		{
			pointphi_2d *temp = tree_point_phi_2d_tail;
			while(temp!=NULL)
			{
				pointphi_2d *temp2 = temp->before;
				delete temp;
				temp = temp2;
			}
		}
		if(tree_grid_2d!=NULL)
		{
			delete_all_quadtree(tree_grid_2d);
		}

		if(mat_savingtree!=NULL)
		{
			int i;
			for(i=0; i<MAXSAVING; i++) delete [] mat_savingtree[i];
			delete [] mat_savingtree;
		}

		//int i;
		
		/*if(level_phi_2d!=NULL)
		{
			for(i=0; i<mesh_x_num+2; i++) delete [] level_phi_2d[i];
			delete [] level_phi_2d;
		}
		
		if(level_burning_rate_2d!=NULL)
		{
			for(i=0; i<mesh_x_num; i++) delete [] level_burning_rate_2d[i];
			delete [] level_burning_rate_2d;
		}*/
	}
	if(dimension == 3)
	{

		if(tree_point_phi_3d_start!=NULL)
		{
			pointphi_3d *temp = tree_point_phi_3d_tail;
			while(temp!=NULL)
			{
				pointphi_3d *temp2 = temp->before;
				delete temp;
				temp = temp2;
			}
		}
		if(tree_grid_3d!=NULL)
		{
			delete_all_octree(tree_grid_3d);
		}

		if(array_savingtree!=NULL) delete [] array_savingtree;

		/*int i; int j;
		
		if(level_phi_3d!=NULL)
		{
			for(i=0; i<mesh_x_num+2; i++)
			{
				for(j=0; j<mesh_y_num+2; j++) delete [] level_phi_3d[i][j];
				delete [] level_phi_3d[i];
			}
			delete [] level_phi_3d;
		}
		
		if(level_burning_rate_3d!=NULL)
		{
			for(i=0; i<mesh_x_num; i++)
			{
				for(j=0; j<mesh_y_num; j++) delete [] level_burning_rate_3d[i][j];
				delete [] level_burning_rate_3d[i];
			}
			delete [] level_burning_rate_3d;
		}*/
	}
}




/*double level_set::tree_hg_old(struct pointphi_2d *point, bool pm, double si1, double si2, double sj1, double sj2)
{
	double *r = quadtree_coeffs(tree_grid_2d,point);
	double *r1 = NULL;
	double *r2 = NULL;
	double *r3 = NULL;
	double *r4 = NULL;

	if(point->left!=NULL) r1 = quadtree_coeffs(tree_grid_2d,point->left);
	if(point->right!=NULL) r4 = quadtree_coeffs(tree_grid_2d,point->right);
	if(point->top!=NULL) r3 = quadtree_coeffs(tree_grid_2d,point->top);
	if(point->bottom!=NULL) r2 = quadtree_coeffs(tree_grid_2d,point->bottom);

	double a;
	double b;
	double c;
	double d;

	if(pm==true)
	{
		a = quadtree_dxplusphi(r,r4,si1);
		b = quadtree_dxminusphi(r,r1,si2);
		c = quadtree_dyplusphi(r,r3,sj1);
		d = quadtree_dyminusphi(r,r2,sj2);
	}
	else
	{
		a = quadtree_dxminusphi(r,r4,si2);
		b = quadtree_dxplusphi(r,r1,si1);
		c = quadtree_dyminusphi(r,r3,sj2);
		d = quadtree_dyplusphi(r,r2,sj1);
	}

	if(r!=NULL) delete r;
	if(r1!=NULL) delete r1;
	if(r2!=NULL) delete r2;
	if(r3!=NULL) delete r3;
	if(r4!=NULL) delete r4;

	return sqrt(max(pow(abs(min(a,0)),2),pow(abs(max(b,0)),2)) + max(pow(abs(min(c,0)),2),pow(abs(max(d,0)),2)));
}

void level_set::iter_reinitial_scheme_2d_old()
{
	pointphi_2d *temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->tempphi3 = temp->phi;
		double *r = quadtree_coeffs(tree_grid_2d,temp);		// phi0, phi1, s1, phi4, s4, phi3, s3, phi2, s2
		if(r[0]*r[3]<0)
		{
			double *r2 = quadtree_coeffs(tree_grid_2d, temp->right);

			double c2 = 1./2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
			double c1 = (r[3] - r[0])/r[4];
			double c0 = (r[3] + r[0])/2. - c2 * pow(r[4],2)/4.;

			if(abs(c2)<MINERROR) temp->si1 = r[4]/2. + (-c0/c1);
			else if(r[0]<0) temp->si1 = r[4]/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->si1 = r[4]/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);

			if(r2!=NULL) delete r2;
		}
		else temp->si1 = 0;
		if(r[0]*r[1]<0)
		{
			double *r2 = quadtree_coeffs(tree_grid_2d, temp->left);

			double c2 = 1./2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
			double c1 = (r[0] - r[1])/r[2];
			double c0 = (r[0] + r[1])/2. - c2 * pow(r[2],2)/4.;

			if(abs(c2)<MINERROR) temp->si2 = r[2]/2. + (-c0/c1);
			else if(r[0]>0) temp->si2 = r[2]/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->si2 = r[2]/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);

			if(r2!=NULL) delete r2;
		}
		else temp->si2 = 0;
		if(r[0]*r[5]<0)
		{
			double *r2 = quadtree_coeffs(tree_grid_2d, temp->top);

			double c2 = 1./2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
			double c1 = (r[5] - r[0])/r[6];
			double c0 = (r[5] + r[0])/2. - c2 * pow(r[6],2)/4.;

			if(abs(c2)<MINERROR) temp->sj1 = r[6]/2. + (-c0/c1);
			else if(r[0]<0) temp->sj1 = r[6]/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->sj1 = r[6]/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);

			if(r2!=NULL) delete r2;
		}
		else temp->sj1 = 0;
		if(r[0]*r[7]<0)
		{
			double *r2 = quadtree_coeffs(tree_grid_2d, temp->bottom);

			double c2 = 1./2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
			double c1 = (r[0] - r[7])/r[8];
			double c0 = (r[7] + r[0])/2. - c2 * pow(r[8],2)/4.;

			if(abs(c2)<MINERROR) temp->sj2 = r[8]/2. + (-c0/c1);
			else if(r[0]>0) temp->sj2 = r[8]/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->sj2 = r[8]/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);

			if(r2!=NULL) delete r2;
		}
		else temp->sj2 = 0;
		if(r!=NULL) delete r;
		temp = temp->next;
	}
	for(int i=0; i<35; i++) reinitial_scheme_2d();
}

void level_set::reinitial_scheme_2d_old()
{
	pointphi_2d *temp = tree_point_phi_2d_start;
	
	while(temp!=NULL)
	{
		double *r = quadtree_coeffs(tree_grid_2d,temp);
		double s1 = 0;
		double s4 = 0;
		double s3 = 0;
		double s2 = 0;
		if(r!=NULL)
		{
			s1 = r[2];
			s4 = r[4];
			s3 = r[6];
			s2 = r[8];

			if(s1==0) s1 = MAXSIZE*10;
			if(s2==0) s2 = MAXSIZE*10;
			if(s3==0) s3 = MAXSIZE*10;
			if(s4==0) s4 = MAXSIZE*10;

			delete r;
		}
		double a = temp->si1; if(a==0) a = MAXSIZE*10;
		double b = temp->si2; if(b==0) b = MAXSIZE*10;
		double c = temp->sj1; if(c==0) c = MAXSIZE*10;
		double d = temp->sj2; if(d==0) d = MAXSIZE*10;
		double si = min(min(min(a,b),c),d);
		double tau = min(min(min(min(s1,s2),s3),s4),si)/2.;
		
		bool pm;
		if(sign(temp->tempphi3)<=0) pm = false;
		else pm = true;

		temp->tempphi = temp->phi - sign(temp->tempphi3) * tau * (tree_hg(temp, pm, temp->si1, temp->si2, temp->sj1, temp->sj2)-1);
		temp = temp->next;
	}
	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->tempphi2 = temp->phi;
		temp->phi = temp->tempphi;
		temp = temp->next;
	}
	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		double *r = quadtree_coeffs(tree_grid_2d,temp);
		double s1 = 0;
		double s4 = 0;
		double s3 = 0;
		double s2 = 0;
		if(r!=NULL)
		{
			s1 = r[2];
			s4 = r[4];
			s3 = r[6];
			s2 = r[8];

			if(s1==0) s1 = MAXSIZE*10;
			if(s2==0) s2 = MAXSIZE*10;
			if(s3==0) s3 = MAXSIZE*10;
			if(s4==0) s4 = MAXSIZE*10;

			delete r;
		}
		double a = temp->si1; if(a==0) a = MAXSIZE*10;
		double b = temp->si2; if(b==0) b = MAXSIZE*10;
		double c = temp->sj1; if(c==0) c = MAXSIZE*10;
		double d = temp->sj2; if(d==0) d = MAXSIZE*10;
		double si = min(min(min(a,b),c),d);
		double tau = min(min(min(min(s1,s2),s3),s4),si)/2.;

		bool pm;
		if(sign(temp->tempphi3)<=0) pm = false;
		else pm = true;

		temp->tempphi = temp->tempphi2 - sign(temp->tempphi3) * tau * (tree_hg(temp, pm, temp->si1, temp->si2, temp->sj1, temp->sj2)-1);
		temp = temp->next;
	}
	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->phi = (temp->tempphi + temp->tempphi2)/2.;
		temp = temp->next;
	}
}

double level_set::g_HJ(int i, int j, int k)
{
	double g;
	if(dimension == 2)
	{
		g = pow(min((level_phi_2d[i][j]-level_phi_2d[i-1][j])/(double)mesh_size, 0.),2);
		g += pow(max((level_phi_2d[i+1][j]-level_phi_2d[i][j])/(double)mesh_size, 0.),2);
		g += pow(min((level_phi_2d[i][j]-level_phi_2d[i][j-1])/(double)mesh_size, 0.),2);
		g += pow(max((level_phi_2d[i][j+1]-level_phi_2d[i][j])/(double)mesh_size, 0.),2);
		g = sqrt(g);
	}
	if(dimension == 3)
	{
		g = pow(min((level_phi_3d[i][j][k]-level_phi_3d[i-1][j][k])/(double)mesh_size, 0.),2);
		g += pow(max((level_phi_3d[i+1][j][k]-level_phi_3d[i][j][k])/(double)mesh_size, 0.),2);
		g += pow(min((level_phi_3d[i][j][k]-level_phi_3d[i][j-1][k])/(double)mesh_size, 0.),2);
		g += pow(max((level_phi_3d[i][j+1][k]-level_phi_3d[i][j][k])/(double)mesh_size, 0.),2);
		g += pow(min((level_phi_3d[i][j][k]-level_phi_3d[i][j][k-1])/(double)mesh_size, 0.),2);
		g += pow(max((level_phi_3d[i][j][k+1]-level_phi_3d[i][j][k])/(double)mesh_size, 0.),2);
		g = sqrt(g);
	}

	return g;
}

bool level_set::in_bounding(int i, int j, int k)
{
	int i0;
	if(dimension==2)
	{
		for(i0=0; i0<level_surface->edge_num; i0++)
		{
			double *a = level_surface->point[level_surface->edge[i0][0]];
			double *b = level_surface->point[level_surface->edge[i0][1]];
			double x1 = min(a[0],b[0]);
			double x2 = max(a[0],b[0]);
			double y1 = min(a[1],b[1]);
			double y2 = max(a[1],b[1]);

			int n1 = (int) floor(x1/mesh_size) - padding_num;
			int n2 = (int) ceil(x2/mesh_size) + padding_num;
			int m1 = (int) floor(y1/mesh_size) - padding_num;
			int m2 = (int) ceil(y2/mesh_size) + padding_num;

			if(n1<-1) n1 = -1;
			if(n2>mesh_x_num) n2 = mesh_x_num;
			if(m1<-1) m1 = -1;
			if(m2>mesh_y_num) m2 = mesh_y_num;

			if(n1+1<i && i<n2+1 && m1+1<j && j<m2+1) return true;
		}
		return false;
	}

	if(dimension==3)
	{
		for(i0=0; i0<level_surface->surf_num; i0++)
		{
			double *a = level_surface->point[level_surface->surf[i0][0]];
			double *b = level_surface->point[level_surface->surf[i0][1]];
			double *c = level_surface->point[level_surface->surf[i0][2]];
			double x1 = min(min(a[0],b[0]),c[0]);
			double x2 = max(max(a[0],b[0]),c[0]);
			double y1 = min(min(a[1],b[1]),c[1]);
			double y2 = max(max(a[1],b[1]),c[1]);
			double z1 = min(min(a[2],b[2]),c[2]);
			double z2 = max(max(a[2],b[2]),c[2]);

			int n1 = (int) floor(x1/mesh_size) - padding_num;
			int n2 = (int) ceil(x2/mesh_size) + padding_num;
			int m1 = (int) floor(y1/mesh_size) - padding_num;
			int m2 = (int) ceil(y2/mesh_size) + padding_num;
			int l1 = (int) floor(z1/mesh_size) - padding_num;
			int l2 = (int) ceil(z2/mesh_size) + padding_num;

			if(n1<-1) n1 = -1;
			if(n2>mesh_x_num) n2 = mesh_x_num;
			if(m1<-1) m1 = -1;
			if(m2>mesh_y_num) m2 = mesh_y_num;
			if(l1<-1) l1 = -1;
			if(l2>mesh_z_num) l2 = mesh_z_num;

			if(n1+1<i && i<n2+1 && m1+1<j && j<m2+1 && l1+1<k && k<l2+1) return true;
		}
		return false;
	}

	return 0;
}

void level_set::bounding_construct_phi(int i0)
{
	int i;
	int j;
	int k;

	if(dimension == 2)
	{
		double *a = level_surface->point[level_surface->edge[i0][0]];
		double *b = level_surface->point[level_surface->edge[i0][1]];
		double x1 = min(a[0],b[0]);
		double x2 = max(a[0],b[0]);
		double y1 = min(a[1],b[1]);
		double y2 = max(a[1],b[1]);

		int n1 = (int) floor(x1/mesh_size) - padding_num;
		int n2 = (int) ceil(x2/mesh_size) + padding_num;
		int m1 = (int) floor(y1/mesh_size) - padding_num;
		int m2 = (int) ceil(y2/mesh_size) + padding_num;

		if(n1<-1) n1 = -1;
		if(n2>mesh_x_num) n2 = mesh_x_num;
		if(m1<-1) m1 = -1;
		if(m2>mesh_y_num) m2 = mesh_y_num;

		for(i=n1; i<=n2; i++)
		{
			for(j=m1; j<=m2; j++)
			{
				if(level_phi_2d[i+1][j+1]!=0)
				{
					double vec[] = {mesh_size * i, mesh_size * j};
					double dist = unsigned_distance_face_point(level_surface, i0, vec);
					if(abs(level_phi_2d[i+1][j+1]) > MAXSIZE * 5)
					{
						if(dist==0) level_phi_2d[i+1][j+1] = 0;
						else level_phi_2d[i+1][j+1] = line_surface_intersecting(level_surface,vec) * dist;
					}
					else if(dist < abs(level_phi_2d[i+1][j+1]))
					{
						if(level_phi_2d[i+1][j+1] > 0) level_phi_2d[i+1][j+1] = dist;
						else level_phi_2d[i+1][j+1] = -dist;
					}
				}
			}
		}
	}

	if(dimension == 3)
	{
		double *a = level_surface->point[level_surface->surf[i0][0]];
		double *b = level_surface->point[level_surface->surf[i0][1]];
		double *c = level_surface->point[level_surface->surf[i0][2]];
		double x1 = min(min(a[0],b[0]),c[0]);
		double x2 = max(max(a[0],b[0]),c[0]);
		double y1 = min(min(a[1],b[1]),c[1]);
		double y2 = max(max(a[1],b[1]),c[1]);
		double z1 = min(min(a[2],b[2]),c[2]);
		double z2 = max(max(a[2],b[2]),c[2]);

		int n1 = (int) floor(x1/mesh_size) - padding_num;
		int n2 = (int) ceil(x2/mesh_size) + padding_num;
		int m1 = (int) floor(y1/mesh_size) - padding_num;
		int m2 = (int) ceil(y2/mesh_size) + padding_num;
		int l1 = (int) floor(z1/mesh_size) - padding_num;
		int l2 = (int) ceil(z2/mesh_size) + padding_num;

		if(n1<-1) n1 = -1;
		if(n2>mesh_x_num) n2 = mesh_x_num;
		if(m1<-1) m1 = -1;
		if(m2>mesh_y_num) m2 = mesh_y_num;
		if(l1<-1) l1 = -1;
		if(l2>mesh_z_num) l2 = mesh_z_num;

		for(i=n1; i<=n2; i++)
		{
			for(j=m1; j<=m2; j++)
			{
				for(k=l1; k<l2; k++)
				{
					if(level_phi_3d[i+1][j+1][k+1]!=0)
					{
						double vec[] = {mesh_size * i, mesh_size * j, mesh_size * k};
						double dist = unsigned_distance_face_point(level_surface, i0, vec);
						if(abs(level_phi_3d[i+1][j+1][k+1]) > MAXSIZE * 5)
						{
							if(dist==0) level_phi_3d[i+1][j+1][k+1] = 0;
							else level_phi_3d[i+1][j+1][k+1] = line_surface_intersecting(level_surface,vec) * dist;
						}
						else if(dist < abs(level_phi_3d[i+1][j+1][k+1]))
						{
							if(level_phi_3d[i+1][j+1][k+1] > 0) level_phi_3d[i+1][j+1][k+1] = dist;
							else level_phi_3d[i+1][j+1][k+1] = -dist;
						}
					}
				}
			}
		}
	}
}

void level_set::construct_phi()
{
	int i;
	int j;
	int k;

	if(dimension == 2)
	{
		if(level_phi_2d==NULL)
		{
			level_phi_2d = new double *[mesh_x_num+2];
			for(i=0; i<mesh_x_num+2; i++)
			{
				level_phi_2d[i] = new double[mesh_y_num+2];
				for(j=0; j<mesh_y_num+2; j++) level_phi_2d[i][j] = MAXSIZE * 10;
			}
		}

		for(i=0; i<level_surface->edge_num; i++) bounding_construct_phi(i);
	}
	if(dimension == 3)
	{
		if(level_phi_3d==NULL)
		{
			level_phi_3d = new double **[mesh_x_num+2];
			for(i=0; i<mesh_x_num+2; i++)
			{
				level_phi_3d[i] = new double *[mesh_y_num+2];
				for(j=0; j<mesh_y_num+2; j++)
				{
					level_phi_3d[i][j] = new double[mesh_z_num+2];
					for(k=0; k<mesh_z_num+2; k++) level_phi_3d[i][j][k] = MAXSIZE * 10;
				}
			}
		}

		for(i=0; i<level_surface->surf_num; i++) bounding_construct_phi(i);
	}
}

void level_set::construct_phi_old()
{
	int i;
	int j;
	int k;

	if(dimension == 2)
	{
		if(level_phi_2d==NULL)
		{
			level_phi_2d = new double *[mesh_x_num+2];
			for(i=0; i<mesh_x_num+2; i++) level_phi_2d[i] = new double[mesh_y_num+2];
		}

		for(i=0; i<mesh_x_num+2; i++)
		{
			for(j=0; j<mesh_y_num+2; j++)
			{
				level_phi_2d[i][j] = distance_surf_point(level_surface, mesh_size * (i-1), mesh_size * (j-1), 0);
			}
		}
	}
	if(dimension == 3)
	{
		if(level_phi_3d==NULL)
		{
			level_phi_3d = new double **[mesh_x_num+2];
			for(i=0; i<mesh_x_num+2; i++)
			{
				level_phi_3d[i] = new double *[mesh_y_num+2];
				for(j=0; j<mesh_y_num+2; j++) level_phi_3d[i][j] = new double[mesh_z_num+2];
			}
		}

		for(i=0; i<mesh_x_num+2; i++)
		{
			for(j=0; j<mesh_y_num+2; j++)
			{
				for(k=0; k<mesh_z_num+2; k++)
				{
					level_phi_3d[i][j][k] = distance_surf_point(level_surface, mesh_size * (i-1), mesh_size * (j-1), mesh_size * (k-1));
				}
			}
		}
	}
}

void level_set::extrapolate_burning_rate()
{
	if(level_surface->burning_rate!=NULL)
	{
		int i;
		int j;
		int k;

		if(dimension == 2)
		{
			if(level_burning_rate_2d == NULL)
			{
				level_burning_rate_2d = new double *[mesh_x_num];
				for(i=0; i<mesh_x_num; i++) level_burning_rate_2d[i] = new double[mesh_y_num];
			}

			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					level_burning_rate_2d[i-1][j-1] = burning_surf_point(level_surface, mesh_size * (i-1), mesh_size * (j-1), 0);
				}
			}
		}
		if(dimension == 3)
		{
			if(level_burning_rate_3d == NULL)
			{
				level_burning_rate_3d = new double **[mesh_x_num];
				for(i=0; i<mesh_x_num; i++)
				{
					level_burning_rate_3d[i] = new double *[mesh_y_num];
					for(j=0; j<mesh_y_num; j++) level_burning_rate_3d[i][j] = new double[mesh_z_num];
				}
			}

			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					for(k=1; k<mesh_z_num+1; k++)
					{
						level_burning_rate_3d[i-1][j-1][k-1] = burning_surf_point(level_surface, mesh_size * (i-1), mesh_size * (j-1), mesh_size * (k-1));
					}
				}
			}
		}
	}

}

void level_set::propagate_surface()
{
	int i; int j; int k;
	
	if(dimension == 2)
	{
		if(accuracy_order == 1)
		{
			double ** temp_phi_2d = new double *[mesh_x_num];
			for(i=0; i<mesh_x_num;i++) temp_phi_2d[i] = new double[mesh_y_num];
			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					if(in_bounding(i,j,0))
					{
					if(level_surface->burning_rate==NULL) temp_phi_2d[i-1][j-1] = level_phi_2d[i][j] + time_step*g_HJ(i,j,0);
					else temp_phi_2d[i-1][j-1] = level_phi_2d[i][j] - level_burning_rate_2d[i-1][j-1]*time_step*g_HJ(i,j,0);
					}
				}
			}

			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					level_phi_2d[i][j] = temp_phi_2d[i-1][j-1];
				}
			}

			for(i=0; i<mesh_x_num;i++) delete [] temp_phi_2d[i];
			delete [] temp_phi_2d;
		}

		propagate_num++;
		if(propagate_num > padding_num/2)
		{
			reconstruct_surface();
			construct_phi();
			propagate_num = 0;
		}
	}
	if(dimension == 3)
	{
		if(accuracy_order == 1)
		{
			double *** temp_phi_3d = new double **[mesh_x_num];
			for(i=0; i<mesh_x_num; i++)
			{
				temp_phi_3d[i] = new double *[mesh_y_num];
				for(j=0; j<mesh_y_num; j++) temp_phi_3d[i][j] = new double[mesh_z_num];
			}
			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					for(k=1; k<mesh_z_num+1; k++)
					{
						if(in_bounding(i,j,k))
						{
						if(level_surface->burning_rate==NULL) temp_phi_3d[i-1][j-1][k-1] = level_phi_3d[i][j][k] + time_step*g_HJ(i,j,k);
						else temp_phi_3d[i-1][j-1][k-1] = level_phi_3d[i][j][k] - level_burning_rate_3d[i-1][j-1][k-1]*time_step*g_HJ(i,j,k);
						}
					}
				}
			}

			for(i=1; i<mesh_x_num+1; i++)
			{
				for(j=1; j<mesh_y_num+1; j++)
				{
					for(k=1; k<mesh_z_num+1; k++)
					{
						level_phi_3d[i][j][k] = temp_phi_3d[i-1][j-1][k-1];
					}
				}
			}

			for(i=0; i<mesh_x_num;i++)
			{
				for(j=0; j<mesh_y_num; j++) delete [] temp_phi_3d[i][j];
				delete [] temp_phi_3d[i];
			}
			delete [] temp_phi_3d;
		}

		propagate_num++;
		if(propagate_num >= padding_num/2)
		{
			reconstruct_surface();
			construct_phi();
			propagate_num = 0;
		}
	}

	// higher accuracy is for next year
}

void level_set::reconstruct_surface()
{
	int i; int j;
	int k1; int k2;
	if(dimension==2)
	{
		int max_point_num = 5*mesh_x_num;
		int max_edge_num = 5*mesh_x_num;
		struct surface tempsurf(2, max_point_num, max_edge_num, false);
		int point_num = 0;
		int edge_num = 0;
		for(i=1; i<mesh_x_num; i++)
		{
			for(j=1; j<mesh_y_num; j++)
			{
				if(in_bounding(i,j,0))
				{
				if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]>=0);			//1
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]<=0);		//2
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]<=0)		//3
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]>=0)		//4
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]<=0)		//5
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]>=0)		//6
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]<=0)		//7
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]>=0)		//8
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
					
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]>=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]<=0)		//9
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]>=0)		//10
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]<=0)		//11
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
					
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]>=0)		//12
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]>=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]<=0)		//13
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)-1.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]>=0)		//14
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]>=0 && level_phi_2d[i+1][j+1]<=0)		//15
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1.), mesh_size *(double(j)-1./2.), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				else if(level_phi_2d[i][j]<=0 && level_phi_2d[i+1][j]<=0 && level_phi_2d[i][j+1]<=0 && level_phi_2d[i+1][j+1]>=0)		//16
				{
					// (mesh_size * (i-1/2), mesh_size * (j-1)), (mesh_size * (i-1), mesh_size * (j-1/2)), (mesh_size * i, mesh_size * (j-1/2)), (mesh_size * (i-1/2), mesh_size * j)
					k1 = add_point(&tempsurf, point_num, mesh_size *(double(i)-1./2.), mesh_size *(double(j)), 0);
					if(k1==point_num) point_num++;
					k2 = add_point(&tempsurf, point_num, mesh_size *(double(i)), mesh_size *(double(j)-1./2.), 0);
					if(k2==point_num) point_num++;
					tempsurf.edge[edge_num][0] = k1;
					tempsurf.edge[edge_num][1] = k2;
					edge_num++;
				}
				}
			}
		}
		
		for(i=0; i<level_surface->point_num; i++) delete [] level_surface->point[i];
		delete [] level_surface->point;
		for(i=0; i<level_surface->edge_num; i++) delete [] level_surface->edge[i];
		delete [] level_surface->edge;
		if(level_surface->burning_rate!=NULL) delete [] level_surface->burning_rate;

		level_surface->point_num = point_num;
		level_surface->point = new double *[point_num];
		for(i=0; i<point_num; i++) level_surface->point[i] = new double[dimension];
		level_surface->edge_num = edge_num;
		level_surface->edge = new int *[edge_num];
		for(i=0; i<edge_num; i++) level_surface->edge[i] = new int[2];
		level_surface->burning_rate = NULL;

		for(i=0; i<point_num; i++)
		{
			level_surface->point[i][0] = tempsurf.point[i][0];
			level_surface->point[i][1] = tempsurf.point[i][1];
		}

		for(i=0; i<edge_num; i++)
		{
			level_surface->edge[i][0] = tempsurf.edge[i][0];
			level_surface->edge[i][1] = tempsurf.edge[i][1];
		}
	}
}*/