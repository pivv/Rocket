#include "levelsets.h"

void level_set::modify_coeffs(double x1, double x2, double y1, double y2)
{
	struct pointphi_2d *temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		if((temp->x==x1 && temp->y>=y1 && temp->y<=y2) || (temp->x==x2 && temp->y>=y1 && temp->y<=y2)
			|| (temp->y==y1 && temp->x>=x1 && temp->x<=x2) || (temp->y==y2 && temp->x>=x1 && temp->x<=x2))
			quadtree_coeffs_saving(tree_grid_2d, temp);
		temp = temp->next;
	}
}

bool level_set::lipschitz_cond(struct quadtree *current, bool b_initial)
{
	double x = current->phi_lefttop->x;
	double y = current->phi_lefttop->y;
	double l = current->length;

	double d1;
	double d2;
	double d3;
	double d4;

	if(b_initial==true)
	{
		/*d1 = abs(distance_circle_point(x,y,0.15,0.5,0.5));
		d2 = abs(distance_circle_point(x,y-l,0.15,0.5,0.5));
		d3 = abs(distance_circle_point(x+l,y,0.15,0.5,0.5));
		d4 = abs(distance_circle_point(x+l,y-l,0.15,0.5,0.5));*/

		d1 = abs(distance_surf_point(level_surface,x,y,0));
		d2 = abs(distance_surf_point(level_surface,x,y-l,0));
		d3 = abs(distance_surf_point(level_surface,x+l,y,0));
		d4 = abs(distance_surf_point(level_surface,x+l,y-l,0));
	}
	else
	{
		d1 = abs(current->phi_lefttop->phi);
		d2 = abs(current->phi_leftbottom->phi);
		d3 = abs(current->phi_righttop->phi);
		d4 = abs(current->phi_rightbottom->phi);
	}

	double d = min(min(min(d1,d2),d3),d4);
	
	if(d < 2.*l*sqrt(2.)) return true;
	else return false;
}

bool level_set::lipschitz_cond(struct octree *current, bool b_initial)
{
	double gamma = 3.;

	double x = current->phi_leftbottomback->x;
	double y = current->phi_leftbottomback->y;
	double z = current->phi_leftbottomback->z;
	double l = current->length;

	double d;

	if(b_initial==true)
	{
		d = abs(distance_surf_point(level_surface,x,y,z));
		d = min(d, abs(distance_surf_point(level_surface,x+l,y,z)));
		d = min(d, abs(distance_surf_point(level_surface,x,y+l,z)));
		d = min(d, abs(distance_surf_point(level_surface,x+l,y+l,z)));
		d = min(d, abs(distance_surf_point(level_surface,x,y,z+l)));
		d = min(d, abs(distance_surf_point(level_surface,x+l,y,z+l)));
		d = min(d, abs(distance_surf_point(level_surface,x,y+l,z+l)));
		d = min(d, abs(distance_surf_point(level_surface,x+l,y+l,z+l)));
	}
	else
	{
		d = abs(current->phi_leftbottomback->phi);
		d = min(d, abs(current->phi_rightbottomback->phi));
		d = min(d, abs(current->phi_lefttopback->phi));
		d = min(d, abs(current->phi_righttopback->phi));
		d = min(d, abs(current->phi_leftbottomfront->phi));
		d = min(d, abs(current->phi_rightbottomfront->phi));
		d = min(d, abs(current->phi_lefttopfront->phi));
		d = min(d, abs(current->phi_righttopfront->phi));
	}
	
	if(d < gamma * l*sqrt(3.)) return true;
	else return false;
}

void level_set::generate_tree(struct quadtree *current, bool b_initial)
{
	if(current == NULL);
	else if(current->tree_lefttop!=NULL)
	{
		if(b_initial==false && lipschitz_cond(current, b_initial)==false)  deletetree(current);
			//&& current->tree_lefttop->tree_lefttop==NULL && current->tree_leftbottom->tree_lefttop==NULL && current->tree_righttop->tree_lefttop==NULL && current->tree_rightbottom->tree_lefttop==NULL) deletetree(current);
		else
		{
			generate_tree(current->tree_lefttop, b_initial);
			generate_tree(current->tree_leftbottom, b_initial);
			generate_tree(current->tree_righttop, b_initial);
			generate_tree(current->tree_rightbottom, b_initial);
		}
	}
	else if(current->depth < max_tree_depth && lipschitz_cond(current, b_initial)==true)
	{
		addtree(current, b_initial);

		generate_tree(current->tree_lefttop, b_initial);
		generate_tree(current->tree_leftbottom, b_initial);
		generate_tree(current->tree_righttop, b_initial);
		generate_tree(current->tree_rightbottom, b_initial);
	}
}

void level_set::generate_tree(struct octree *current, bool b_initial)
{
	if(current == NULL);
	else if(current->tree_leftbottomback!=NULL)
	{
		if(b_initial==false && lipschitz_cond(current, b_initial)==false)  deletetree(current);
		else
		{
			generate_tree(current->tree_leftbottomback, b_initial);
			generate_tree(current->tree_lefttopback, b_initial);
			generate_tree(current->tree_rightbottomback, b_initial);
			generate_tree(current->tree_righttopback, b_initial);
			generate_tree(current->tree_leftbottomfront, b_initial);
			generate_tree(current->tree_lefttopfront, b_initial);
			generate_tree(current->tree_rightbottomfront, b_initial);
			generate_tree(current->tree_righttopfront, b_initial);
		}
	}
	else if(current->depth < max_tree_depth && lipschitz_cond(current, b_initial)==true)
	{
		addtree(current, b_initial);
		
		generate_tree(current->tree_leftbottomback, b_initial);
		generate_tree(current->tree_lefttopback, b_initial);
		generate_tree(current->tree_rightbottomback, b_initial);
		generate_tree(current->tree_righttopback, b_initial);
		generate_tree(current->tree_leftbottomfront, b_initial);
		generate_tree(current->tree_lefttopfront, b_initial);
		generate_tree(current->tree_rightbottomfront, b_initial);
		generate_tree(current->tree_righttopfront, b_initial);
	}
}

void level_set::deletetree(struct quadtree *current)
{
	if(current->tree_lefttop->tree_lefttop!=NULL) deletetree(current->tree_lefttop);
	if(current->tree_leftbottom->tree_lefttop!=NULL) deletetree(current->tree_leftbottom);
	if(current->tree_righttop->tree_lefttop!=NULL) deletetree(current->tree_righttop);
	if(current->tree_rightbottom->tree_lefttop!=NULL) deletetree(current->tree_rightbottom);

	double x = current->phi_lefttop->x;
	double y = current->phi_lefttop->y;
	double l = current->length;

	deletelink(current->tree_lefttop->phi_rightbottom);
	current->tree_lefttop->phi_leftbottom->right = NULL;
	current->tree_righttop->phi_rightbottom->left = NULL;
	current->tree_lefttop->phi_righttop->bottom = NULL;
	current->tree_leftbottom->phi_rightbottom->top = NULL;

	if(current->tree_lefttop->phi_leftbottom->left==NULL) deletelink(current->tree_lefttop->phi_leftbottom);
	if(current->tree_righttop->phi_rightbottom->right==NULL) deletelink(current->tree_righttop->phi_rightbottom);
	if(current->tree_lefttop->phi_righttop->top==NULL) deletelink(current->tree_lefttop->phi_righttop);
	if(current->tree_leftbottom->phi_rightbottom->bottom==NULL) deletelink(current->tree_leftbottom->phi_rightbottom);
	
	delete current->tree_lefttop;
	delete current->tree_leftbottom;
	delete current->tree_righttop;
	delete current->tree_rightbottom;

	current->tree_lefttop = NULL;
	current->tree_leftbottom = NULL;
	current->tree_righttop = NULL;
	current->tree_rightbottom = NULL;

	//modify_coeffs(x,x+l,y-l,y);

	tree_num-=4;
	tree_end_num-=3;
}

void level_set::deletetree(struct octree *current)
{
	if(current->tree_leftbottomback->tree_leftbottomback!=NULL) deletetree(current->tree_leftbottomback);
	if(current->tree_lefttopback->tree_leftbottomback!=NULL) deletetree(current->tree_lefttopback);
	if(current->tree_rightbottomback->tree_leftbottomback!=NULL) deletetree(current->tree_rightbottomback);
	if(current->tree_righttopback->tree_leftbottomback!=NULL) deletetree(current->tree_righttopback);
	if(current->tree_leftbottomfront->tree_leftbottomback!=NULL) deletetree(current->tree_leftbottomfront);
	if(current->tree_lefttopfront->tree_leftbottomback!=NULL) deletetree(current->tree_lefttopfront);
	if(current->tree_rightbottomfront->tree_leftbottomback!=NULL) deletetree(current->tree_rightbottomfront);
	if(current->tree_righttopfront->tree_leftbottomback!=NULL) deletetree(current->tree_righttopfront);

	deletelink(current->tree_leftbottomback->phi_righttopfront, false);

	if(current->tree_leftbottomback->phi_lefttopfront->left == NULL) deletelink(current->tree_leftbottomback->phi_lefttopfront, false);
	if(current->tree_leftbottomback->phi_rightbottomfront->bottom == NULL) deletelink(current->tree_leftbottomback->phi_rightbottomfront, false);
	if(current->tree_leftbottomback->phi_righttopback->back == NULL) deletelink(current->tree_leftbottomback->phi_righttopback, false);
	if(current->tree_righttopfront->phi_rightbottomback->right == NULL) deletelink(current->tree_righttopfront->phi_rightbottomback, false);
	if(current->tree_righttopfront->phi_lefttopback->top == NULL) deletelink(current->tree_righttopfront->phi_lefttopback, false);
	if(current->tree_righttopfront->phi_leftbottomfront->front == NULL) deletelink(current->tree_righttopfront->phi_leftbottomfront, false);
	
	//phi_leftbottom
	if(current->tree_leftbottomback->phi_leftbottomfront->left == NULL
		&& current->tree_leftbottomback->phi_leftbottomfront->bottom == NULL) deletelink(current->tree_leftbottomback->phi_leftbottomfront, true);
	//phi_bottomback
	if(current->tree_leftbottomback->phi_rightbottomback->bottom == NULL
		&& current->tree_leftbottomback->phi_rightbottomback->back == NULL) deletelink(current->tree_leftbottomback->phi_rightbottomback, true);
	//phi_leftback
	if(current->tree_leftbottomback->phi_lefttopback->left == NULL
		&& current->tree_leftbottomback->phi_lefttopback->back == NULL) deletelink(current->tree_leftbottomback->phi_lefttopback, true);
	
	//phi_righttop
	if(current->tree_righttopback->phi_righttopfront->right == NULL
		&& current->tree_righttopback->phi_righttopfront->top == NULL) deletelink(current->tree_righttopback->phi_righttopfront, true);
	//phi_topback
	if(current->tree_righttopback->phi_lefttopback->top == NULL
		&& current->tree_righttopback->phi_lefttopback->back == NULL) deletelink(current->tree_righttopback->phi_lefttopback, true);
	//phi_rightback
	if(current->tree_righttopback->phi_rightbottomback->right == NULL
		&& current->tree_righttopback->phi_rightbottomback->back == NULL) deletelink(current->tree_righttopback->phi_rightbottomback, true);
	
	//phi_lefttop
	if(current->tree_lefttopfront->phi_lefttopback->left == NULL
		&& current->tree_lefttopfront->phi_lefttopback->top == NULL) deletelink(current->tree_lefttopfront->phi_lefttopback, true);
	//phi_topfront
	if(current->tree_lefttopfront->phi_righttopfront->top == NULL
		&& current->tree_lefttopfront->phi_righttopfront->front == NULL) deletelink(current->tree_lefttopfront->phi_righttopfront, true);
	//phi_leftfront
	if(current->tree_lefttopfront->phi_leftbottomfront->left == NULL
		&& current->tree_lefttopfront->phi_leftbottomfront->front == NULL) deletelink(current->tree_lefttopfront->phi_leftbottomfront, true);
	
	//phi_rightbottom
	if(current->tree_rightbottomfront->phi_rightbottomback->right == NULL
		&& current->tree_rightbottomfront->phi_rightbottomback->bottom == NULL) deletelink(current->tree_rightbottomfront->phi_rightbottomback, true);
	//phi_bottomfront
	if(current->tree_rightbottomfront->phi_leftbottomfront->bottom == NULL
		&& current->tree_rightbottomfront->phi_leftbottomfront->front == NULL) deletelink(current->tree_rightbottomfront->phi_leftbottomfront, true);
	//phi_rightfront
	if(current->tree_rightbottomfront->phi_righttopfront->right == NULL
		&& current->tree_rightbottomfront->phi_righttopfront->front == NULL) deletelink(current->tree_rightbottomfront->phi_righttopfront, true);

	delete current->tree_leftbottomback;
	delete current->tree_lefttopback;
	delete current->tree_rightbottomback;
	delete current->tree_righttopback;
	delete current->tree_leftbottomfront;
	delete current->tree_lefttopfront;
	delete current->tree_rightbottomfront;
	delete current->tree_righttopfront;
	
	current->tree_leftbottomback = NULL;
	current->tree_lefttopback = NULL;
	current->tree_rightbottomback = NULL;
	current->tree_righttopback = NULL;
	current->tree_leftbottomfront = NULL;
	current->tree_lefttopfront = NULL;
	current->tree_rightbottomfront = NULL;
	current->tree_righttopfront = NULL;

	tree_num-=8;
	tree_end_num-=7;
}

void level_set::addtree(struct quadtree *current, bool b_initial)
{
	double x = current->phi_lefttop->x;
	double y = current->phi_lefttop->y;
	double l = current->length;

	/*struct pointphi_2d *phi_left = find_vertex(tree_grid_2d, x, y-l/2.);
	struct pointphi_2d *phi_right = find_vertex(tree_grid_2d, x+l, y-l/2.);
	struct pointphi_2d *phi_top = find_vertex(tree_grid_2d, x+l/2., y);
	struct pointphi_2d *phi_bottom = find_vertex(tree_grid_2d, x+l/2., y-l);*/

	struct pointphi_2d *phi_left = find_vertex(current->phi_lefttop, 1, x, y-l/2.);
	struct pointphi_2d *phi_right = find_vertex(current->phi_righttop, 1, x+l, y-l/2.);
	struct pointphi_2d *phi_bottom = find_vertex(current->phi_rightbottom, 0, x+l/2., y-l);
	struct pointphi_2d *phi_top = find_vertex(current->phi_righttop, 0, x+l/2., y);
	struct pointphi_2d *phi_center;

	if(b_initial==true)
	{
		/*if(phi_left==NULL) phi_left = addlink(x, y-l/2., distance_circle_point(x,y-l/2.,0.15,0.5,0.5), NULL, NULL, current->phi_lefttop, current->phi_leftbottom);
		if(phi_right==NULL) phi_right = addlink(x+l, y-l/2., distance_circle_point(x+l,y-l/2.,0.15,0.5,0.5), NULL, NULL, current->phi_righttop, current->phi_rightbottom);
		if(phi_top==NULL) phi_top = addlink(x+l/2., y, distance_circle_point(x+l/2.,y,0.15,0.5,0.5), current->phi_lefttop, current->phi_righttop, NULL, NULL);
		if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y-l, distance_circle_point(x+l/2.,y-l,0.15,0.5,0.5), current->phi_leftbottom, current->phi_rightbottom, NULL, NULL);
		phi_center = addlink(x+l/2., y-l/2., distance_circle_point(x+l/2.,y-l/2.,0.15,0.5,0.5), phi_left, phi_right, phi_top, phi_bottom);*/

		if(phi_left==NULL) phi_left = addlink(x, y-l/2., distance_surf_point(level_surface,x,y-l/2.,0), NULL, NULL, current->phi_lefttop, current->phi_leftbottom);
		if(phi_right==NULL) phi_right = addlink(x+l, y-l/2., distance_surf_point(level_surface,x+l,y-l/2.,0), NULL, NULL, current->phi_righttop, current->phi_rightbottom);
		if(phi_top==NULL) phi_top = addlink(x+l/2., y, distance_surf_point(level_surface,x+l/2.,y,0), current->phi_lefttop, current->phi_righttop, NULL, NULL);
		if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y-l, distance_surf_point(level_surface,x+l/2.,y-l,0), current->phi_leftbottom, current->phi_rightbottom, NULL, NULL);
		phi_center = addlink(x+l/2., y-l/2., distance_surf_point(level_surface,x+l/2.,y-l/2.,0), phi_left, phi_right, phi_top, phi_bottom);
	}
	else
	{
		if(phi_left==NULL) phi_left = addlink(x, y-l/2., computephi(current,x,y-l/2.), NULL, NULL, current->phi_lefttop, current->phi_leftbottom);
		if(phi_right==NULL) phi_right = addlink(x+l, y-l/2., computephi(current,x+l,y-l/2.), NULL, NULL, current->phi_righttop, current->phi_rightbottom);
		if(phi_top==NULL) phi_top = addlink(x+l/2., y, computephi(current,x+l/2.,y), current->phi_lefttop, current->phi_righttop, NULL, NULL);
		if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y-l, computephi(current,x+l/2.,y-l), current->phi_leftbottom, current->phi_rightbottom, NULL, NULL);
		phi_center = addlink(x+l/2., y-l/2., computephi(current,x+l/2.,y-l/2.), phi_left, phi_right, phi_top, phi_bottom);
	}

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

	//quadtree_coeffs_saving(tree_grid_2d, phi_center);
	//modify_coeffs(x,x+l,y-l,y);
	
	tree_num+=4;
	tree_end_num+=3;
}

void level_set::addtree(struct octree *current, bool b_initial)
{
	double x = current->phi_leftbottomback->x;
	double y = current->phi_leftbottomback->y;
	double z = current->phi_leftbottomback->z;
	double l = current->length;

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

	/*struct pointphi_3d *phi_left = find_vertex(tree_grid_3d, x, y+l/2., z+l/2.);
	struct pointphi_3d *phi_right = find_vertex(tree_grid_3d, x+l, y+l/2., z+l/2.);
	struct pointphi_3d *phi_bottom = find_vertex(tree_grid_3d, x+l/2., y, z+l/2.);
	struct pointphi_3d *phi_top = find_vertex(tree_grid_3d, x+l/2., y+l, z+l/2.);
	struct pointphi_3d *phi_back = find_vertex(tree_grid_3d, x+l/2., y+l/2., z);
	struct pointphi_3d *phi_front = find_vertex(tree_grid_3d, x+l/2., y+l/2., z+l);

	struct pointphi_3d *phi_center;

	struct pointphi_3d *phi_leftbottom = find_vertex(tree_grid_3d, x, y, z+l/2.);
	struct pointphi_3d *phi_rightbottom = find_vertex(tree_grid_3d, x+l, y, z+l/2.);
	struct pointphi_3d *phi_lefttop = find_vertex(tree_grid_3d, x, y+l, z+l/2.);
	struct pointphi_3d *phi_righttop = find_vertex(tree_grid_3d, x+l, y+l, z+l/2.);

	struct pointphi_3d *phi_leftback = find_vertex(tree_grid_3d, x, y+l/2., z);
	struct pointphi_3d *phi_leftfront = find_vertex(tree_grid_3d, x, y+l/2., z+l);
	struct pointphi_3d *phi_rightback = find_vertex(tree_grid_3d, x+l, y+l/2., z);
	struct pointphi_3d *phi_rightfront = find_vertex(tree_grid_3d, x+l, y+l/2., z+l);

	struct pointphi_3d *phi_bottomback = find_vertex(tree_grid_3d, x+l/2., y, z);
	struct pointphi_3d *phi_bottomfront = find_vertex(tree_grid_3d, x+l/2., y, z+l);
	struct pointphi_3d *phi_topback = find_vertex(tree_grid_3d, x+l/2., y+l, z);
	struct pointphi_3d *phi_topfront = find_vertex(tree_grid_3d, x+l/2., y+l, z+l);*/

	if(b_initial==true)
	{
		if(phi_leftbottom==NULL) phi_leftbottom = addlink(x, y, z+l/2., distance_surf_point(level_surface, x, y, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_leftbottomback, current->phi_leftbottomfront);
		if(phi_rightbottom==NULL) phi_rightbottom = addlink(x+l, y, z+l/2., distance_surf_point(level_surface, x+l, y, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_rightbottomback, current->phi_rightbottomfront);
		if(phi_lefttop==NULL) phi_lefttop = addlink(x, y+l, z+l/2., distance_surf_point(level_surface, x, y+l, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_lefttopback, current->phi_lefttopfront);
		if(phi_righttop==NULL) phi_righttop = addlink(x+l, y+l, z+l/2., distance_surf_point(level_surface, x+l, y+l, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_righttopback, current->phi_righttopfront);
		
		if(phi_leftback==NULL) phi_leftback = addlink(x, y+l/2., z, distance_surf_point(level_surface, x, y+l/2., z),
			NULL, NULL, current->phi_leftbottomback, current->phi_lefttopback, NULL, NULL);
		if(phi_leftfront==NULL) phi_leftfront = addlink(x, y+l/2., z+l, distance_surf_point(level_surface, x, y+l/2., z+l),
			NULL, NULL, current->phi_leftbottomfront, current->phi_lefttopfront, NULL, NULL);
		if(phi_rightback==NULL) phi_rightback = addlink(x+l, y+l/2., z, distance_surf_point(level_surface, x+l, y+l/2., z),
			NULL, NULL, current->phi_rightbottomback, current->phi_righttopback, NULL, NULL);
		if(phi_rightfront==NULL) phi_rightfront = addlink(x+l, y+l/2., z+l, distance_surf_point(level_surface, x+l, y+l/2., z+l),
			NULL, NULL, current->phi_rightbottomfront, current->phi_righttopfront, NULL, NULL);
		
		if(phi_bottomback==NULL) phi_bottomback = addlink(x+l/2., y, z, distance_surf_point(level_surface, x+l/2., y, z),
			current->phi_leftbottomback, current->phi_rightbottomback, NULL, NULL, NULL, NULL);
		if(phi_bottomfront==NULL) phi_bottomfront = addlink(x+l/2., y, z+l, distance_surf_point(level_surface, x+l/2., y, z+l),
			current->phi_leftbottomfront, current->phi_rightbottomfront, NULL, NULL, NULL, NULL);
		if(phi_topback==NULL) phi_topback = addlink(x+l/2., y+l, z, distance_surf_point(level_surface, x+l/2., y+l, z),
			current->phi_lefttopback, current->phi_righttopback, NULL, NULL, NULL, NULL);
		if(phi_topfront==NULL) phi_topfront = addlink(x+l/2., y+l, z+l, distance_surf_point(level_surface, x+l/2., y+l, z+l),
			current->phi_lefttopfront, current->phi_righttopfront, NULL, NULL, NULL, NULL);

		if(phi_left==NULL) phi_left = addlink(x, y+l/2., z+l/2., distance_surf_point(level_surface, x, y+l/2., z+l/2.), NULL, NULL, phi_leftbottom, phi_lefttop, phi_leftback, phi_leftfront);
		if(phi_right==NULL) phi_right = addlink(x+l, y+l/2., z+l/2., distance_surf_point(level_surface, x+l, y+l/2., z+l/2.), NULL, NULL, phi_rightbottom, phi_righttop, phi_rightback, phi_rightfront);
		if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y, z+l/2., distance_surf_point(level_surface, x+l/2., y, z+l/2.), phi_leftbottom, phi_rightbottom, NULL, NULL, phi_bottomback,phi_bottomfront);
		if(phi_top==NULL) phi_top = addlink(x+l/2., y+l, z+l/2., distance_surf_point(level_surface, x+l/2., y+l, z+l/2.), phi_lefttop, phi_righttop, NULL, NULL, phi_topback, phi_topfront);
		if(phi_back==NULL) phi_back = addlink(x+l/2., y+l/2., z, distance_surf_point(level_surface, x+l/2., y+l/2., z), phi_leftback, phi_rightback, phi_bottomback, phi_topback, NULL, NULL);
		if(phi_front==NULL) phi_front = addlink(x+l/2., y+l/2., z+l, distance_surf_point(level_surface, x+l/2., y+l/2., z+l), phi_leftfront, phi_rightfront, phi_bottomfront, phi_topfront, NULL, NULL);

		phi_center = addlink(x+l/2., y+l/2., z+l/2., distance_surf_point(level_surface, x+l/2., y+l/2., z+l/2.), phi_left, phi_right, phi_bottom, phi_top, phi_back, phi_front);
	}
	else
	{
		if(phi_leftbottom==NULL) phi_leftbottom = addlink(x, y, z+l/2., computephi(current, x, y, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_leftbottomback, current->phi_leftbottomfront);
		if(phi_rightbottom==NULL) phi_rightbottom = addlink(x+l, y, z+l/2., computephi(current, x+l, y, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_rightbottomback, current->phi_rightbottomfront);
		if(phi_lefttop==NULL) phi_lefttop = addlink(x, y+l, z+l/2., computephi(current, x, y+l, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_lefttopback, current->phi_lefttopfront);
		if(phi_righttop==NULL) phi_righttop = addlink(x+l, y+l, z+l/2., computephi(current, x+l, y+l, z+l/2.),
			NULL, NULL, NULL, NULL, current->phi_righttopback, current->phi_righttopfront);
		
		if(phi_leftback==NULL) phi_leftback = addlink(x, y+l/2., z, computephi(current, x, y+l/2., z),
			NULL, NULL, current->phi_leftbottomback, current->phi_lefttopback, NULL, NULL);
		if(phi_leftfront==NULL) phi_leftfront = addlink(x, y+l/2., z+l, computephi(current, x, y+l/2., z+l),
			NULL, NULL, current->phi_leftbottomfront, current->phi_lefttopfront, NULL, NULL);
		if(phi_rightback==NULL) phi_rightback = addlink(x+l, y+l/2., z, computephi(current, x+l, y+l/2., z),
			NULL, NULL, current->phi_rightbottomback, current->phi_righttopback, NULL, NULL);
		if(phi_rightfront==NULL) phi_rightfront = addlink(x+l, y+l/2., z+l, computephi(current, x+l, y+l/2., z+l),
			NULL, NULL, current->phi_rightbottomfront, current->phi_righttopfront, NULL, NULL);
		
		if(phi_bottomback==NULL) phi_bottomback = addlink(x+l/2., y, z, computephi(current, x+l/2., y, z),
			current->phi_leftbottomback, current->phi_rightbottomback, NULL, NULL, NULL, NULL);
		if(phi_bottomfront==NULL) phi_bottomfront = addlink(x+l/2., y, z+l, computephi(current, x+l/2., y, z+l),
			current->phi_leftbottomfront, current->phi_rightbottomfront, NULL, NULL, NULL, NULL);
		if(phi_topback==NULL) phi_topback = addlink(x+l/2., y+l, z, computephi(current, x+l/2., y+l, z),
			current->phi_lefttopback, current->phi_righttopback, NULL, NULL, NULL, NULL);
		if(phi_topfront==NULL) phi_topfront = addlink(x+l/2., y+l, z+l, computephi(current, x+l/2., y+l, z+l),
			current->phi_lefttopfront, current->phi_righttopfront, NULL, NULL, NULL, NULL);

		if(phi_left==NULL) phi_left = addlink(x, y+l/2., z+l/2., computephi(current, x, y+l/2., z+l/2.), NULL, NULL, phi_leftbottom, phi_lefttop, phi_leftback, phi_leftfront);
		if(phi_right==NULL) phi_right = addlink(x+l, y+l/2., z+l/2., computephi(current, x+l, y+l/2., z+l/2.), NULL, NULL, phi_rightbottom, phi_righttop, phi_rightback, phi_rightfront);
		if(phi_bottom==NULL) phi_bottom = addlink(x+l/2., y, z+l/2., computephi(current, x+l/2., y, z+l/2.), phi_leftbottom, phi_rightbottom, NULL, NULL, phi_bottomback,phi_bottomfront);
		if(phi_top==NULL) phi_top = addlink(x+l/2., y+l, z+l/2., computephi(current, x+l/2., y+l, z+l/2.), phi_lefttop, phi_righttop, NULL, NULL, phi_topback, phi_topfront);
		if(phi_back==NULL) phi_back = addlink(x+l/2., y+l/2., z, computephi(current, x+l/2., y+l/2., z), phi_leftback, phi_rightback, phi_bottomback, phi_topback, NULL, NULL);
		if(phi_front==NULL) phi_front = addlink(x+l/2., y+l/2., z+l, computephi(current, x+l/2., y+l/2., z+l), phi_leftfront, phi_rightfront, phi_bottomfront, phi_topfront, NULL, NULL);

		phi_center = addlink(x+l/2., y+l/2., z+l/2., computephi(current, x+l/2., y+l/2., z+l/2.), phi_left, phi_right, phi_bottom, phi_top, phi_back, phi_front);
	}


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

	tree_num+=8;
	tree_end_num+=7;
}
