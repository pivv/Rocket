#include "memoryuse.h"
#include "basictools.h"

struct pointphi_2d
{
	double x;
	double y;
	double phi;
	double tempphi;
	double tempphi2;
	double tempphi3;
	double si1;
	double si2;
	double sj1;
	double sj2;

	struct pointphi_2d *left;
	struct pointphi_2d *right;
	struct pointphi_2d *top;
	struct pointphi_2d *bottom;

	double leftphi;
	double rightphi;
	double topphi;
	double bottomphi;
	double lefts;
	double rights;
	double tops;
	double bottoms;

	struct pointphi_2d *next;
	struct pointphi_2d *before;
	pointphi_2d();
};

pointphi_2d::pointphi_2d()
{
	left = NULL;
	right = NULL;
	top = NULL;
	bottom = NULL;
	next = NULL;
	before = NULL;
}

struct pointphi_3d
{
	double x;
	double y;
	double z;
	double phi;
	double tempphi;

	double tempphi2;
	double tempphi3;
	
	/*double si1;
	double si2;
	double sj1;
	double sj2;
	double sk1;
	double sk2;*/

	struct pointphi_3d *left;
	struct pointphi_3d *right;
	struct pointphi_3d *top;
	struct pointphi_3d *bottom;
	struct pointphi_3d *front;
	struct pointphi_3d *back;

	/*double leftphi;
	double rightphi;
	double topphi;
	double bottomphi;
	double frontphi;
	double backphi;
	double lefts;
	double rights;
	double tops;
	double bottoms;
	double fronts;
	double backs;*/

	struct pointphi_3d *next;
	struct pointphi_3d *before;
	pointphi_3d();
};

pointphi_3d::pointphi_3d()
{
	left = NULL;
	right = NULL;
	top = NULL;
	bottom = NULL;
	front = NULL;
	back = NULL;
	next = NULL;
	before = NULL;
}

struct quadtree
{
	int depth;

	double length;

	struct pointphi_2d *phi_lefttop;
	struct pointphi_2d *phi_righttop;
	struct pointphi_2d *phi_leftbottom;
	struct pointphi_2d *phi_rightbottom;

	struct quadtree *tree_lefttop;
	struct quadtree *tree_righttop;
	struct quadtree *tree_leftbottom;
	struct quadtree *tree_rightbottom;

	struct quadtree *parent;

	quadtree();
};

quadtree::quadtree()
{
	phi_lefttop = NULL;
	phi_righttop = NULL;
	phi_leftbottom = NULL;
	phi_rightbottom = NULL;
	tree_lefttop = NULL;
	tree_righttop = NULL;
	tree_leftbottom = NULL;
	tree_rightbottom = NULL;
	parent = NULL;
}

struct octree
{
	int depth;

	double length;

	struct pointphi_3d *phi_lefttopfront;
	struct pointphi_3d *phi_righttopfront;
	struct pointphi_3d *phi_leftbottomfront;
	struct pointphi_3d *phi_rightbottomfront;
	struct pointphi_3d *phi_lefttopback;
	struct pointphi_3d *phi_righttopback;
	struct pointphi_3d *phi_leftbottomback;
	struct pointphi_3d *phi_rightbottomback;

	struct octree *tree_lefttopfront;
	struct octree *tree_righttopfront;
	struct octree *tree_leftbottomfront;
	struct octree *tree_rightbottomfront;
	struct octree *tree_lefttopback;
	struct octree *tree_righttopback;
	struct octree *tree_leftbottomback;
	struct octree *tree_rightbottomback;

	struct octree *parent;

	octree();

	octree(struct octree *tree_parent, int tree_depth, double tree_length, struct pointphi_3d *tree_phi_leftbottomback,
		struct pointphi_3d *tree_phi_lefttopback, struct pointphi_3d *tree_phi_rightbottomback, struct pointphi_3d *tree_phi_righttopback,
		struct pointphi_3d *tree_phi_leftbottomfront, struct pointphi_3d *tree_phi_lefttopfront, struct pointphi_3d *tree_phi_rightbottomfront,
		struct pointphi_3d *tree_phi_righttopfront);
};

octree::octree()
{
	phi_leftbottomback = NULL;
	phi_rightbottomback = NULL;
	phi_lefttopback = NULL;
	phi_righttopback = NULL;
	phi_leftbottomfront = NULL;
	phi_rightbottomfront = NULL;
	phi_lefttopfront = NULL;
	phi_righttopfront = NULL;

	tree_leftbottomback = NULL;
	tree_rightbottomback = NULL;
	tree_lefttopback = NULL;
	tree_righttopback = NULL;
	tree_leftbottomfront = NULL;
	tree_rightbottomfront = NULL;
	tree_lefttopfront = NULL;
	tree_righttopfront = NULL;

	parent = NULL;
}

octree::octree(struct octree *tree_parent, int tree_depth, double tree_length, struct pointphi_3d *tree_phi_leftbottomback,
		struct pointphi_3d *tree_phi_lefttopback, struct pointphi_3d *tree_phi_rightbottomback, struct pointphi_3d *tree_phi_righttopback,
		struct pointphi_3d *tree_phi_leftbottomfront, struct pointphi_3d *tree_phi_lefttopfront, struct pointphi_3d *tree_phi_rightbottomfront,
		struct pointphi_3d *tree_phi_righttopfront)
{
	depth = tree_depth;
	length = tree_length;

	phi_leftbottomback = tree_phi_leftbottomback;
	phi_rightbottomback = tree_phi_rightbottomback;
	phi_lefttopback = tree_phi_lefttopback;
	phi_righttopback = tree_phi_righttopback;
	phi_leftbottomfront = tree_phi_leftbottomfront;
	phi_rightbottomfront = tree_phi_rightbottomfront;
	phi_lefttopfront = tree_phi_lefttopfront;
	phi_righttopfront = tree_phi_righttopfront;
	
	tree_leftbottomback = NULL;
	tree_rightbottomback = NULL;
	tree_lefttopback = NULL;
	tree_righttopback = NULL;
	tree_leftbottomfront = NULL;
	tree_rightbottomfront = NULL;
	tree_lefttopfront = NULL;
	tree_righttopfront = NULL;

	parent = tree_parent;
}

void delete_all_quadtree(struct quadtree *current)
{
	if(current!=NULL)
	{
		delete_all_quadtree(current->tree_lefttop);
		delete_all_quadtree(current->tree_leftbottom);
		delete_all_quadtree(current->tree_righttop);
		delete_all_quadtree(current->tree_rightbottom);

		delete current;
	}
}

void delete_all_octree(struct octree *current)
{
	if(current!=NULL)
	{
		delete_all_octree(current->tree_leftbottomback);
		delete_all_octree(current->tree_lefttopback);
		delete_all_octree(current->tree_rightbottomback);
		delete_all_octree(current->tree_righttopback);
		delete_all_octree(current->tree_leftbottomfront);
		delete_all_octree(current->tree_lefttopfront);
		delete_all_octree(current->tree_rightbottomfront);
		delete_all_octree(current->tree_righttopfront);

		delete current;
	}
}

struct quadtree *find_left_top_face(struct quadtree *current, double x, double y)
{
	struct quadtree *temp = current;
	while(temp->tree_lefttop!=NULL)
	{
		if(x<=temp->phi_lefttop->x + temp->length/2. + MINERROR)
		{
			if(y>=temp->phi_lefttop->y - temp->length/2. - MINERROR) temp = temp->tree_lefttop;
			else temp = temp->tree_leftbottom;
		}
		else
		{
			if(y>=temp->phi_lefttop->y - temp->length/2. - MINERROR) temp = temp->tree_righttop;
			else temp = temp->tree_rightbottom;
		}
	}

	return temp;
}

struct quadtree *find_right_bottom_face(struct quadtree *current, double x, double y)
{
	struct quadtree *temp = current;
	while(temp->tree_lefttop!=NULL)
	{
		if(x<temp->phi_lefttop->x + temp->length/2. - MINERROR)
		{
			if(y>temp->phi_lefttop->y - temp->length/2. + MINERROR) temp = temp->tree_lefttop;
			else temp = temp->tree_leftbottom;
		}
		else
		{
			if(y>temp->phi_lefttop->y - temp->length/2. + MINERROR) temp = temp->tree_righttop;
			else temp = temp->tree_rightbottom;
		}
	}

	return temp;
}


struct octree *find_left_bottom_back_face(struct octree *current, double x, double y, double z)
{
	struct octree *temp = current;
	while(temp->tree_leftbottomback!=NULL)
	{
		if(x<=temp->phi_leftbottomback->x + temp->length/2. + MINERROR)
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. + MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_leftbottomback;
				else temp = temp->tree_leftbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_lefttopback;
				else temp = temp->tree_lefttopfront;
			}
		}
		else
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. + MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_rightbottomback;
				else temp = temp->tree_rightbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_righttopback;
				else temp = temp->tree_righttopfront;
			}
		}
	}

	return temp;
}


struct octree *find_right_top_back_face(struct octree *current, double x, double y, double z)
{
	struct octree *temp = current;
	while(temp->tree_leftbottomback!=NULL)
	{
		if(x<=temp->phi_leftbottomback->x + temp->length/2. - MINERROR)
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. - MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_leftbottomback;
				else temp = temp->tree_leftbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_lefttopback;
				else temp = temp->tree_lefttopfront;
			}
		}
		else
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. - MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_rightbottomback;
				else temp = temp->tree_rightbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. + MINERROR) temp = temp->tree_righttopback;
				else temp = temp->tree_righttopfront;
			}
		}
	}

	return temp;
}


struct octree *find_left_bottom_front_face(struct octree *current, double x, double y, double z)
{
	struct octree *temp = current;
	while(temp->tree_leftbottomback!=NULL)
	{
		if(x<=temp->phi_leftbottomback->x + temp->length/2. + MINERROR)
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. + MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_leftbottomback;
				else temp = temp->tree_leftbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_lefttopback;
				else temp = temp->tree_lefttopfront;
			}
		}
		else
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. + MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_rightbottomback;
				else temp = temp->tree_rightbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_righttopback;
				else temp = temp->tree_righttopfront;
			}
		}
	}

	return temp;
}


struct octree *find_right_top_front_face(struct octree *current, double x, double y, double z)
{
	struct octree *temp = current;
	while(temp->tree_leftbottomback!=NULL)
	{
		if(x<=temp->phi_leftbottomback->x + temp->length/2. - MINERROR)
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. - MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_leftbottomback;
				else temp = temp->tree_leftbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_lefttopback;
				else temp = temp->tree_lefttopfront;
			}
		}
		else
		{
			if(y<=temp->phi_leftbottomback->y + temp->length/2. - MINERROR)
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_rightbottomback;
				else temp = temp->tree_rightbottomfront;
			}
			else
			{
				if(z<=temp->phi_leftbottomback->z + temp->length/2. - MINERROR) temp = temp->tree_righttopback;
				else temp = temp->tree_righttopfront;
			}
		}
	}

	return temp;
}

struct pointphi_2d *find_vertex(struct quadtree *current, double x, double y)
{
	struct quadtree *temp = find_left_top_face(current,x,y);

	if(abs(x-temp->phi_lefttop->x) < MINERROR)
	{
		if(abs(y-temp->phi_lefttop->y) < MINERROR) return temp->phi_lefttop;
		else if(abs(y-temp->phi_leftbottom->y) < MINERROR) return temp->phi_leftbottom;
	}
	else if(abs(x-temp->phi_righttop->x) < MINERROR)
	{
		if(abs(y-temp->phi_righttop->y) < MINERROR) return temp->phi_righttop;
		else if(abs(y-temp->phi_rightbottom->y) < MINERROR) return temp->phi_rightbottom;
	}

	temp = find_right_bottom_face(current,x,y);

	if(abs(x-temp->phi_lefttop->x) < MINERROR)
	{
		if(abs(y-temp->phi_lefttop->y) < MINERROR) return temp->phi_lefttop;
		else if(abs(y-temp->phi_leftbottom->y) < MINERROR) return temp->phi_leftbottom;
	}
	else if(abs(x-temp->phi_righttop->x) < MINERROR)
	{
		if(abs(y-temp->phi_righttop->y) < MINERROR) return temp->phi_righttop;
		else if(abs(y-temp->phi_rightbottom->y) < MINERROR) return temp->phi_rightbottom;
	}

	return NULL;
}

struct pointphi_3d *find_vertex(struct octree *current, double x, double y, double z)
{
	struct octree *temp = find_left_bottom_back_face(current,x,y,z);

	if(abs(x-temp->phi_leftbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_leftbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_leftbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_lefttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_lefttopfront;
		}
	}
	else if(abs(x-temp->phi_rightbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_rightbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_rightbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_righttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_righttopfront;
		}
	}

	temp = find_right_top_back_face(current,x,y,z);

	if(abs(x-temp->phi_leftbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_leftbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_leftbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_lefttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_lefttopfront;
		}
	}
	else if(abs(x-temp->phi_rightbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_rightbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_rightbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_righttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_righttopfront;
		}
	}
	
	temp = find_left_bottom_front_face(current,x,y,z);

	if(abs(x-temp->phi_leftbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_leftbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_leftbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_lefttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_lefttopfront;
		}
	}
	else if(abs(x-temp->phi_rightbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_rightbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_rightbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_righttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_righttopfront;
		}
	}

	temp = find_right_top_front_face(current,x,y,z);

	if(abs(x-temp->phi_leftbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_leftbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_leftbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_lefttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_lefttopfront;
		}
	}
	else if(abs(x-temp->phi_rightbottomback->x) < MINERROR)
	{
		if(abs(y-temp->phi_leftbottomback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_rightbottomback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_rightbottomfront;
		}
		if(abs(y-temp->phi_lefttopback->y) < MINERROR)
		{
			if(abs(z-temp->phi_leftbottomback->z) < MINERROR) return temp->phi_righttopback;
			else if(abs(z-temp->phi_leftbottomfront->z) < MINERROR) return temp->phi_righttopfront;
		}
	}


	return NULL;
}