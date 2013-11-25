#include "levelsets_8test.h"

struct pointphi_2d *level_set::addlink(double x, double y, double phi, struct pointphi_2d *left, struct pointphi_2d *right, struct pointphi_2d *top, struct pointphi_2d *bottom, int type)
{
	if(type==0)
	{
		struct pointphi_2d *newpoint = new struct pointphi_2d();
		if(tree_rcase_point_phi_2d_tail!=NULL) tree_rcase_point_phi_2d_tail->next = newpoint;
		newpoint->before = tree_rcase_point_phi_2d_tail;
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
		newpoint->b_rate = 1.;

		if(tree_rcase_point_phi_2d_start==NULL) tree_rcase_point_phi_2d_start = newpoint;

		tree_rcase_point_phi_2d_tail = newpoint;

		link_num++;

		return newpoint;
	}
	if(type==1)
	{
		struct pointphi_2d *newpoint = new struct pointphi_2d();
		if(tree_propel_point_phi_2d_tail!=NULL) tree_propel_point_phi_2d_tail->next = newpoint;
		newpoint->before = tree_propel_point_phi_2d_tail;
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
		newpoint->b_rate = 1.;

		if(tree_propel_point_phi_2d_start==NULL) tree_propel_point_phi_2d_start = newpoint;

		tree_propel_point_phi_2d_tail = newpoint;

		link_num++;

		return newpoint;
	}
	return NULL;
}

struct pointphi_3d *level_set::addlink(double x, double y, double z, double phi, struct pointphi_3d *left, struct pointphi_3d *right, struct pointphi_3d *bottom, struct pointphi_3d *top, struct pointphi_3d *back, struct pointphi_3d *front, int type)
{
	if(type==0)
	{
		struct pointphi_3d *newpoint = new struct pointphi_3d();
		if(tree_rcase_point_phi_3d_tail!=NULL) tree_rcase_point_phi_3d_tail->next = newpoint;
		newpoint->before = tree_rcase_point_phi_3d_tail;
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
		newpoint->b_rate = 1.;

		if(tree_rcase_point_phi_3d_start==NULL) tree_rcase_point_phi_3d_start = newpoint;

		tree_rcase_point_phi_3d_tail = newpoint;

		link_num++;

		return newpoint;
	}
	if(type==1)
	{
		struct pointphi_3d *newpoint = new struct pointphi_3d();
		if(tree_propel_point_phi_3d_tail!=NULL) tree_propel_point_phi_3d_tail->next = newpoint;
		newpoint->before = tree_propel_point_phi_3d_tail;
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
		newpoint->b_rate = 1.;

		if(tree_propel_point_phi_3d_start==NULL) tree_propel_point_phi_3d_start = newpoint;

		tree_propel_point_phi_3d_tail = newpoint;

		link_num++;

		return newpoint;
	}
	return NULL;
}

void level_set::initialize_tree(int type)
{
	if(dimension==2)
	{
		struct surface *surf;
		if(type==0)
		{
			tree_rcase_grid_2d = new quadtree();
			tree_rcase_grid_2d->depth = 0;
			tree_rcase_grid_2d->length = (double)MAXSIZE;
			surf = rcase_surface;
		}
		if(type==1)
		{
			tree_propel_grid_2d = new quadtree();
			tree_propel_grid_2d->depth = 0;
			tree_propel_grid_2d->length = (double)MAXSIZE;
			surf = propel_surface;
		}
 
		struct pointphi_2d *initleftbottom = addlink(0.,0.,distance_surf_point(surf,0.,0.,0), NULL, NULL, NULL, NULL, type);
		struct pointphi_2d *initlefttop = addlink(0.,(double)MAXSIZE,distance_surf_point(surf,0.,(double)MAXSIZE,0), NULL, NULL, NULL, NULL, type);
		struct pointphi_2d *initrightbottom = addlink((double)MAXSIZE,0.,distance_surf_point(surf,(double)MAXSIZE,0.,0), NULL, NULL, NULL, NULL, type);
		struct pointphi_2d *initrighttop = addlink((double)MAXSIZE,(double)MAXSIZE,distance_surf_point(surf,(double)MAXSIZE,(double)MAXSIZE,0), NULL, NULL, NULL, NULL, type);

		initleftbottom->right = initrightbottom;
		initleftbottom->top = initlefttop;

		initlefttop->right = initrighttop;
		initlefttop->bottom = initleftbottom;

		initrightbottom->left = initleftbottom;
		initrightbottom->top = initrighttop;

		initrighttop->left = initlefttop;
		initrighttop->bottom = initleftbottom;

		if(type==0)
		{
			tree_rcase_grid_2d->phi_leftbottom = initleftbottom;
			tree_rcase_grid_2d->phi_lefttop = initlefttop;
			tree_rcase_grid_2d->phi_rightbottom = initrightbottom;
			tree_rcase_grid_2d->phi_righttop = initrighttop;
		}
		if(type==1)
		{
			tree_propel_grid_2d->phi_leftbottom = initleftbottom;
			tree_propel_grid_2d->phi_lefttop = initlefttop;
			tree_propel_grid_2d->phi_rightbottom = initrightbottom;
			tree_propel_grid_2d->phi_righttop = initrighttop;
		}
	}
	if(dimension==3)
	{
		struct surface *surf;

		if(type==0)
		{
			tree_rcase_grid_3d = new octree();
			tree_rcase_grid_3d->depth = 0;
			tree_rcase_grid_3d->length = (double)MAXSIZE;
			surf = rcase_surface;
		}
		if(type==1)
		{
			tree_propel_grid_3d = new octree();
			tree_propel_grid_3d->depth = 0;
			tree_propel_grid_3d->length = (double)MAXSIZE;
			surf = propel_surface;
		}

		struct pointphi_3d *initleftbottomback = addlink(0.,0.,0.,distance_surf_point(surf,0.,0.,0.), NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initlefttopback = addlink(0.,(double)MAXSIZE,0.,distance_surf_point(surf,0.,(double)MAXSIZE,0.),NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initrightbottomback = addlink((double)MAXSIZE,0.,0.,distance_surf_point(surf,(double)MAXSIZE,0.,0.),NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initrighttopback = addlink((double)MAXSIZE,(double)MAXSIZE,0.,distance_surf_point(surf,(double)MAXSIZE,(double)MAXSIZE,0.),NULL, NULL, NULL, NULL, NULL, NULL, type);
		
		struct pointphi_3d *initleftbottomfront = addlink(0.,0.,(double)MAXSIZE,distance_surf_point(surf,0.,0.,(double)MAXSIZE), NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initlefttopfront = addlink(0.,(double)MAXSIZE,(double)MAXSIZE,distance_surf_point(surf,0.,(double)MAXSIZE,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initrightbottomfront = addlink((double)MAXSIZE,0.,(double)MAXSIZE,distance_surf_point(surf,(double)MAXSIZE,0.,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL, type);
		struct pointphi_3d *initrighttopfront = addlink((double)MAXSIZE,(double)MAXSIZE,(double)MAXSIZE,distance_surf_point(surf,(double)MAXSIZE,(double)MAXSIZE,(double)MAXSIZE),NULL, NULL, NULL, NULL, NULL, NULL, type);


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
		
		if(type==0)
		{
			tree_rcase_grid_3d->phi_leftbottomback = initleftbottomback;
			tree_rcase_grid_3d->phi_lefttopback = initlefttopback;
			tree_rcase_grid_3d->phi_rightbottomback = initrightbottomback;
			tree_rcase_grid_3d->phi_righttopback = initrighttopback;

			tree_rcase_grid_3d->phi_leftbottomfront = initleftbottomfront;
			tree_rcase_grid_3d->phi_lefttopfront = initlefttopfront;
			tree_rcase_grid_3d->phi_rightbottomfront = initrightbottomfront;
			tree_rcase_grid_3d->phi_righttopfront = initrighttopfront;
		}
		if(type==1)
		{
			tree_propel_grid_3d->phi_leftbottomback = initleftbottomback;
			tree_propel_grid_3d->phi_lefttopback = initlefttopback;
			tree_propel_grid_3d->phi_rightbottomback = initrightbottomback;
			tree_propel_grid_3d->phi_righttopback = initrighttopback;

			tree_propel_grid_3d->phi_leftbottomfront = initleftbottomfront;
			tree_propel_grid_3d->phi_lefttopfront = initlefttopfront;
			tree_propel_grid_3d->phi_rightbottomfront = initrightbottomfront;
			tree_propel_grid_3d->phi_righttopfront = initrighttopfront;
		}
	}
}


bool level_set::lipschitz_cond(struct quadtree *current, int type)
{
	double x = current->phi_lefttop->x;
	double y = current->phi_lefttop->y;
	double l = current->length;

	double d1;
	double d2;
	double d3;
	double d4;

	struct surface *surf;
	if(type==0)
	{
		surf = rcase_surface;
	}
	if(type==1)
	{
		surf = propel_surface;
	}

	d1 = abs(distance_surf_point(surf,x,y,0));
	d2 = abs(distance_surf_point(surf,x,y-l,0));
	d3 = abs(distance_surf_point(surf,x+l,y,0));
	d4 = abs(distance_surf_point(surf,x+l,y-l,0));

	double d = min(min(min(d1,d2),d3),d4);
	
	if(d < 1.*l*sqrt(2.)) return true;
	else return false;
}

bool level_set::lipschitz_cond(struct octree *current, int type)
{
	double gamma = 1.;

	double x = current->phi_leftbottomback->x;
	double y = current->phi_leftbottomback->y;
	double z = current->phi_leftbottomback->z;
	double l = current->length;

	double d;

	struct surface *surf;
	if(type==0)
	{
		surf = rcase_surface;
	}
	if(type==1)
	{
		surf = propel_surface;
	}

	d = abs(distance_surf_point(surf,x,y,z));
	d = min(d, abs(distance_surf_point(surf,x+l,y,z)));
	d = min(d, abs(distance_surf_point(surf,x,y+l,z)));
	d = min(d, abs(distance_surf_point(surf,x+l,y+l,z)));
	d = min(d, abs(distance_surf_point(surf,x,y,z+l)));
	d = min(d, abs(distance_surf_point(surf,x+l,y,z+l)));
	d = min(d, abs(distance_surf_point(surf,x,y+l,z+l)));
	d = min(d, abs(distance_surf_point(surf,x+l,y+l,z+l)));
	
	if(d < gamma * l*sqrt(3.)) return true;
	else return false;
}

void level_set::generate_tree(struct quadtree *current, int type)
{
	if(current == NULL);
	else if(current->depth < max_tree_depth && lipschitz_cond(current, type)==true)
	{
		addtree(current, type);

		generate_tree(current->tree_lefttop, type);
		generate_tree(current->tree_leftbottom, type);
		generate_tree(current->tree_righttop, type);
		generate_tree(current->tree_rightbottom, type);
	}
}

void level_set::generate_tree(struct octree *current, int type)
{
	if(current == NULL);
	else if(current->depth < max_tree_depth && lipschitz_cond(current, type)==true)
	{
		addtree(current, type);
		
		generate_tree(current->tree_leftbottomback, type);
		generate_tree(current->tree_lefttopback, type);
		generate_tree(current->tree_rightbottomback, type);
		generate_tree(current->tree_righttopback, type);
		generate_tree(current->tree_leftbottomfront, type);
		generate_tree(current->tree_lefttopfront, type);
		generate_tree(current->tree_rightbottomfront, type);
		generate_tree(current->tree_righttopfront, type);
	}
}


void level_set::addtree(struct quadtree *current, int type)
{
	double x = current->phi_lefttop->x;
	double y = current->phi_lefttop->y;
	double l = current->length;

	struct surface *surf;
	if(type==0)
	{
		surf = rcase_surface;
	}
	if(type==1)
	{
		surf = propel_surface;
	}

	struct pointphi_2d *phi_left = find_vertex(current->phi_lefttop, 1, x, y-l/2.);
	struct pointphi_2d *phi_right = find_vertex(current->phi_righttop, 1, x+l, y-l/2.);
	struct pointphi_2d *phi_bottom = find_vertex(current->phi_rightbottom, 0, x+l/2., y-l);
	struct pointphi_2d *phi_top = find_vertex(current->phi_righttop, 0, x+l/2., y);
	struct pointphi_2d *phi_center;

	if(phi_left==NULL) phi_left = addlink(x, y-l/2., distance_surf_point(surf,x,y-l/2.,0), NULL, NULL, current->phi_lefttop, current->phi_leftbottom, type);
	if(phi_right==NULL) phi_right = addlink(x+l, y-l/2., distance_surf_point(surf,x+l,y-l/2.,0), NULL, NULL, current->phi_righttop, current->phi_rightbottom, type);
	if(phi_top==NULL) phi_top = addlink(x+l/2., y, distance_surf_point(surf,x+l/2.,y,0), current->phi_lefttop, current->phi_righttop, NULL, NULL, type);
	if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y-l, distance_surf_point(surf,x+l/2.,y-l,0), current->phi_leftbottom, current->phi_rightbottom, NULL, NULL, type);
	phi_center = addlink(x+l/2., y-l/2., distance_surf_point(surf,x+l/2.,y-l/2.,0), phi_left, phi_right, phi_top, phi_bottom, type);

	current->tree_lefttop = new struct quadtree;
	current->tree_leftbottom = new struct quadtree;
	current->tree_righttop = new struct quadtree;
	current->tree_rightbottom = new struct quadtree;

	current->tree_lefttop->parent = current;
	current->tree_lefttop->depth = current->depth + 1;
	current->tree_lefttop->length = current->length/2.;
	current->tree_lefttop->phi_lefttop = current->phi_lefttop;
	current->tree_lefttop->phi_leftbottom = phi_left;
	current->tree_lefttop->phi_righttop = phi_top;
	current->tree_lefttop->phi_rightbottom = phi_center;
	
	current->tree_leftbottom->parent = current;
	current->tree_leftbottom->depth = current->depth + 1;
	current->tree_leftbottom->length = current->length/2.;
	current->tree_leftbottom->phi_lefttop = phi_left;
	current->tree_leftbottom->phi_leftbottom = current->phi_leftbottom;
	current->tree_leftbottom->phi_righttop = phi_center;
	current->tree_leftbottom->phi_rightbottom = phi_bottom;
	
	current->tree_righttop->parent = current;
	current->tree_righttop->depth = current->depth + 1;
	current->tree_righttop->length = current->length/2.;
	current->tree_righttop->phi_lefttop = phi_top;
	current->tree_righttop->phi_leftbottom = phi_center;
	current->tree_righttop->phi_righttop = current->phi_righttop;
	current->tree_righttop->phi_rightbottom = phi_right;
	
	current->tree_rightbottom->parent = current;
	current->tree_rightbottom->depth = current->depth + 1;
	current->tree_rightbottom->length = current->length/2.;
	current->tree_rightbottom->phi_lefttop = phi_center;
	current->tree_rightbottom->phi_leftbottom = phi_bottom;
	current->tree_rightbottom->phi_righttop = phi_right;
	current->tree_rightbottom->phi_rightbottom = current->phi_rightbottom;
}

void level_set::addtree(struct octree *current, int type)
{
	double x = current->phi_leftbottomback->x;
	double y = current->phi_leftbottomback->y;
	double z = current->phi_leftbottomback->z;
	double l = current->length;

	struct surface *surf;
	if(type==0)
	{
		surf = rcase_surface;
	}
	if(type==1)
	{
		surf = propel_surface;
	}

	// left bottom back : 0 1 2
	// leftbottom leftback bottomback : 3 4 5

	struct pointphi_3d *phi_left = find_vertex(current->phi_lefttopfront, 5, x, y+l/2., z+l/2.);
	struct pointphi_3d *phi_right = find_vertex(current->phi_righttopfront, 5, x+l, y+l/2., z+l/2.);
	struct pointphi_3d *phi_bottom = find_vertex(current->phi_rightbottomfront, 4, x+l/2., y, z+l/2.);
	struct pointphi_3d *phi_top = find_vertex(current->phi_righttopfront, 4, x+l/2., y+l, z+l/2.);
	struct pointphi_3d *phi_back = find_vertex(current->phi_righttopback, 3, x+l/2., y+l/2., z);
	struct pointphi_3d *phi_front = find_vertex(current->phi_righttopfront, 3, x+l/2., y+l/2., z+l);

	struct pointphi_3d *phi_center;

	struct pointphi_3d *phi_leftbottom = find_vertex(current->phi_leftbottomfront, 2, x, y, z+l/2.);
	struct pointphi_3d *phi_rightbottom = find_vertex(current->phi_rightbottomfront, 2, x+l, y, z+l/2.);
	struct pointphi_3d *phi_lefttop = find_vertex(current->phi_lefttopfront, 2, x, y+l, z+l/2.);
	struct pointphi_3d *phi_righttop = find_vertex(current->phi_righttopfront, 2, x+l, y+l, z+l/2.);

	struct pointphi_3d *phi_leftback = find_vertex(current->phi_lefttopback, 1, x, y+l/2., z);
	struct pointphi_3d *phi_leftfront = find_vertex(current->phi_lefttopfront, 1, x, y+l/2., z+l);
	struct pointphi_3d *phi_rightback = find_vertex(current->phi_righttopback, 1, x+l, y+l/2., z);
	struct pointphi_3d *phi_rightfront = find_vertex(current->phi_righttopfront, 1, x+l, y+l/2., z+l);

	struct pointphi_3d *phi_bottomback = find_vertex(current->phi_rightbottomback, 0, x+l/2., y, z);
	struct pointphi_3d *phi_bottomfront = find_vertex(current->phi_rightbottomfront, 0, x+l/2., y, z+l);
	struct pointphi_3d *phi_topback = find_vertex(current->phi_righttopback, 0, x+l/2., y+l, z);
	struct pointphi_3d *phi_topfront = find_vertex(current->phi_righttopfront, 0, x+l/2., y+l, z+l);

	if(phi_leftbottom==NULL) phi_leftbottom = addlink(x, y, z+l/2., distance_surf_point(surf, x, y, z+l/2.),
		NULL, NULL, NULL, NULL, current->phi_leftbottomback, current->phi_leftbottomfront, type);
	if(phi_rightbottom==NULL) phi_rightbottom = addlink(x+l, y, z+l/2., distance_surf_point(surf, x+l, y, z+l/2.),
		NULL, NULL, NULL, NULL, current->phi_rightbottomback, current->phi_rightbottomfront, type);
	if(phi_lefttop==NULL) phi_lefttop = addlink(x, y+l, z+l/2., distance_surf_point(surf, x, y+l, z+l/2.),
		NULL, NULL, NULL, NULL, current->phi_lefttopback, current->phi_lefttopfront, type);
	if(phi_righttop==NULL) phi_righttop = addlink(x+l, y+l, z+l/2., distance_surf_point(surf, x+l, y+l, z+l/2.),
		NULL, NULL, NULL, NULL, current->phi_righttopback, current->phi_righttopfront, type);
		
	if(phi_leftback==NULL) phi_leftback = addlink(x, y+l/2., z, distance_surf_point(surf, x, y+l/2., z),
		NULL, NULL, current->phi_leftbottomback, current->phi_lefttopback, NULL, NULL, type);
	if(phi_leftfront==NULL) phi_leftfront = addlink(x, y+l/2., z+l, distance_surf_point(surf, x, y+l/2., z+l),
		NULL, NULL, current->phi_leftbottomfront, current->phi_lefttopfront, NULL, NULL, type);
	if(phi_rightback==NULL) phi_rightback = addlink(x+l, y+l/2., z, distance_surf_point(surf, x+l, y+l/2., z),
		NULL, NULL, current->phi_rightbottomback, current->phi_righttopback, NULL, NULL, type);
	if(phi_rightfront==NULL) phi_rightfront = addlink(x+l, y+l/2., z+l, distance_surf_point(surf, x+l, y+l/2., z+l),
		NULL, NULL, current->phi_rightbottomfront, current->phi_righttopfront, NULL, NULL, type);
		
	if(phi_bottomback==NULL) phi_bottomback = addlink(x+l/2., y, z, distance_surf_point(surf, x+l/2., y, z),
		current->phi_leftbottomback, current->phi_rightbottomback, NULL, NULL, NULL, NULL, type);
	if(phi_bottomfront==NULL) phi_bottomfront = addlink(x+l/2., y, z+l, distance_surf_point(surf, x+l/2., y, z+l),
		current->phi_leftbottomfront, current->phi_rightbottomfront, NULL, NULL, NULL, NULL, type);
	if(phi_topback==NULL) phi_topback = addlink(x+l/2., y+l, z, distance_surf_point(surf, x+l/2., y+l, z),
		current->phi_lefttopback, current->phi_righttopback, NULL, NULL, NULL, NULL, type);
	if(phi_topfront==NULL) phi_topfront = addlink(x+l/2., y+l, z+l, distance_surf_point(surf, x+l/2., y+l, z+l),
		current->phi_lefttopfront, current->phi_righttopfront, NULL, NULL, NULL, NULL, type);

	if(phi_left==NULL) phi_left = addlink(x, y+l/2., z+l/2., distance_surf_point(surf, x, y+l/2., z+l/2.), NULL, NULL, phi_leftbottom, phi_lefttop, phi_leftback, phi_leftfront, type);
	if(phi_right==NULL) phi_right = addlink(x+l, y+l/2., z+l/2., distance_surf_point(surf, x+l, y+l/2., z+l/2.), NULL, NULL, phi_rightbottom, phi_righttop, phi_rightback, phi_rightfront, type);
	if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y, z+l/2., distance_surf_point(surf, x+l/2., y, z+l/2.), phi_leftbottom, phi_rightbottom, NULL, NULL, phi_bottomback,phi_bottomfront, type);
	if(phi_top==NULL) phi_top = addlink(x+l/2., y+l, z+l/2., distance_surf_point(surf, x+l/2., y+l, z+l/2.), phi_lefttop, phi_righttop, NULL, NULL, phi_topback, phi_topfront, type);
	if(phi_back==NULL) phi_back = addlink(x+l/2., y+l/2., z, distance_surf_point(surf, x+l/2., y+l/2., z), phi_leftback, phi_rightback, phi_bottomback, phi_topback, NULL, NULL, type);
	if(phi_front==NULL) phi_front = addlink(x+l/2., y+l/2., z+l, distance_surf_point(surf, x+l/2., y+l/2., z+l), phi_leftfront, phi_rightfront, phi_bottomfront, phi_topfront, NULL, NULL, type);

	phi_center = addlink(x+l/2., y+l/2., z+l/2., distance_surf_point(surf, x+l/2., y+l/2., z+l/2.), phi_left, phi_right, phi_bottom, phi_top, phi_back, phi_front, type);

	current->tree_leftbottomback = new struct octree(current, current->depth + 1, current->length/2.,
		current->phi_leftbottomback, phi_leftback, phi_bottomback, phi_back, phi_leftbottom, phi_left, phi_bottom, phi_center);
	
	current->tree_lefttopback = new struct octree(current, current->depth + 1, current->length/2.,
		phi_leftback, current->phi_lefttopback, phi_back, phi_topback, phi_left, phi_lefttop, phi_center, phi_top);

	current->tree_rightbottomback = new struct octree(current, current->depth + 1, current->length/2.,
		phi_bottomback, phi_back, current->phi_rightbottomback, phi_rightback, phi_bottom, phi_center, phi_rightbottom, phi_right);
	
	current->tree_righttopback = new struct octree(current, current->depth + 1, current->length/2.,
		phi_back, phi_topback, phi_rightback, current->phi_righttopback, phi_center, phi_top, phi_right, phi_righttop);
	

	current->tree_leftbottomfront = new struct octree(current, current->depth + 1, current->length/2.,
		phi_leftbottom, phi_left, phi_bottom, phi_center, current->phi_leftbottomfront, phi_leftfront, phi_bottomfront, phi_front);
	
	current->tree_lefttopfront = new struct octree(current, current->depth + 1, current->length/2.,
		phi_left, phi_lefttop, phi_center, phi_top, phi_leftfront, current->phi_lefttopfront, phi_front, phi_topfront);

	current->tree_rightbottomfront = new struct octree(current, current->depth + 1, current->length/2.,
		phi_bottom, phi_center, phi_rightbottom, phi_right, phi_bottomfront, phi_front, current->phi_rightbottomfront, phi_rightfront);
	
	current->tree_righttopfront = new struct octree(current, current->depth + 1, current->length/2.,
		phi_center, phi_top, phi_right, phi_righttop, phi_front, phi_topfront, phi_rightfront, current->phi_righttopfront);
}

struct quadtree *level_set::find_cell_containing_point_type(double x, double y, int type)
{
	struct quadtree *current1;
	struct quadtree *current2;
	if(type==0)
	{
		current1 = find_left_top_face(tree_rcase_grid_2d,x,y);
		current2 = find_right_bottom_face(tree_rcase_grid_2d,x,y);
	}
	if(type==1)
	{
		current1 = find_left_top_face(tree_propel_grid_2d,x,y);
		current2 = find_right_bottom_face(tree_propel_grid_2d,x,y);
	}
	if(current1->length < current2->length)	return current1;
	else return current2;
}

double level_set::level_computephi_type(double x, double y, int type)
{
	return computephi(find_cell_containing_point_type(x,y,type),x,y);
}

struct octree *level_set::find_cell_containing_point_type(double x, double y, double z, int type)
{
	struct octree *current1;
	struct octree *current2;
	struct octree *current3;
	struct octree *current4;

	if(type==0)
	{
		current1 = find_left_bottom_back_face(tree_rcase_grid_3d,x,y,z);
		current2 = find_right_top_back_face(tree_rcase_grid_3d,x,y,z);
		current3 = find_left_bottom_front_face(tree_rcase_grid_3d,x,y,z);
		current4 = find_right_top_front_face(tree_rcase_grid_3d,x,y,z);
	}
	if(type==1)
	{
		current1 = find_left_bottom_back_face(tree_propel_grid_3d,x,y,z);
		current2 = find_right_top_back_face(tree_propel_grid_3d,x,y,z);
		current3 = find_left_bottom_front_face(tree_propel_grid_3d,x,y,z);
		current4 = find_right_top_front_face(tree_propel_grid_3d,x,y,z);
	}

	int mini = 1;
	double minleng = current1->length;
	if(current2->length < minleng)
	{
		mini = 2;
		minleng = current2->length;
	}
	if(current3->length < minleng)
	{
		mini = 3;
		minleng = current3->length;
	}
	if(current4->length < minleng)
	{
		mini = 4;
		minleng = current4->length;
	}

	if(mini==1) return current1;
	else if(mini==2) return current2;
	else if(mini==3) return current3;
	else return current4;
}

double level_set::level_computephi_type(double x, double y, double z, int type)
{
	return computephi(find_cell_containing_point_type(x,y,z,type),x,y,z);
}