#include "levelsets_2extra.h"

struct quadtree *level_set::find_cell_containing_point(double x, double y)
{
	struct quadtree *current1 = find_left_top_face(tree_grid_2d,x,y);
	struct quadtree *current2 = find_right_bottom_face(tree_grid_2d,x,y);
	if(current1->length < current2->length)	return current1;
	else return current2;
}

double level_set::level_computephi(double x, double y)
{
	return computephi(find_cell_containing_point(x,y),x,y);
}

struct octree *level_set::find_cell_containing_point(double x, double y, double z)
{
	struct octree *current1 = find_left_bottom_back_face(tree_grid_3d,x,y,z);
	struct octree *current2 = find_right_top_back_face(tree_grid_3d,x,y,z);
	struct octree *current3 = find_left_bottom_front_face(tree_grid_3d,x,y,z);
	struct octree *current4 = find_right_top_front_face(tree_grid_3d,x,y,z);

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

double level_set::level_computephi(double x, double y, double z)
{
	return computephi(find_cell_containing_point(x,y,z),x,y,z);
}



//here



double external_velocity_x_2d(double x, double y, int num)
{
	return -(y-(double)MAXSIZE/2.);
	if(num < 256) return -pow(sin(PI*x/(double)MAXSIZE),2) * sin(2*PI*y/(double)MAXSIZE);
	else return pow(sin(PI*x/(double)MAXSIZE),2) * sin(2*PI*y/(double)MAXSIZE);
}

double external_velocity_y_2d(double x, double y, int num)
{
	return (x-(double)MAXSIZE/2.);
	if(num < 256) return pow(sin(PI*y/(double)MAXSIZE),2) * sin(2*PI*x/(double)MAXSIZE);
	else return -pow(sin(PI*y/(double)MAXSIZE),2) * sin(2*PI*x/(double)MAXSIZE);
}


double external_velocity_x_3d(double x, double y, double z, int num)
{
	return 0.3;
}

double external_velocity_y_3d(double x, double y, double z, int num)
{
	return 0.3;
}

double external_velocity_z_3d(double x, double y, double z, int num)
{
	return 0.3;
}


double convection_velocity(pointphi_2d *temp, double delta, int num)
{
	double phi_x; double phi_y;
	double u = external_velocity_x_2d(temp->x, temp->y, num);
	double v = external_velocity_y_2d(temp->x, temp->y, num);

	// dxplus dxminus dyplus dyminus dzplus dzminus : 0 1 2 3 4 5
	if(u >= 0) phi_x = quadtree_weno(temp,1,delta);
	else phi_x = quadtree_weno(temp,0,delta);
	if(v >= 0) phi_y = quadtree_weno(temp,3,delta);
	else phi_y = quadtree_weno(temp,2,delta);
	
	/*if(u >= 0) phi_x = quadtree_dxminus(temp);
	else phi_x = quadtree_dxplus(temp);
	if(v >= 0) phi_y = quadtree_dyminus(temp);
	else phi_y = quadtree_dyplus(temp);*/

	return -(u*phi_x + v*phi_y);
}

double convection_velocity(pointphi_3d *temp, double delta, int num)
{
	double phi_x; double phi_y; double phi_z;
	double u = external_velocity_x_3d(temp->x, temp->y, temp->z, num);
	double v = external_velocity_y_3d(temp->x, temp->y, temp->z, num);
	double w = external_velocity_z_3d(temp->x, temp->y, temp->z, num);

	// dxplus dxminus dyplus dyminus dzplus dzminus : 0 1 2 3 4 5
	if(u >= 0) phi_x = octree_weno(temp,1,delta);
	else phi_x = octree_weno(temp,0,delta);
	if(v >= 0) phi_y = octree_weno(temp,3,delta);
	else phi_y = octree_weno(temp,2,delta);
	if(w >= 0) phi_z = octree_weno(temp,5,delta);
	else phi_z = octree_weno(temp,4,delta);

	/*if(u >= 0) phi_x = octree_dxminus(temp);
	else phi_x = octree_dxplus(temp);
	if(v >= 0) phi_y = octree_dyminus(temp);
	else phi_y = octree_dyplus(temp);
	if(w >= 0) phi_z = octree_dzminus(temp);
	else phi_z = octree_dzplus(temp);*/

	return -(u*phi_x + v*phi_y + w*phi_z);
}

void level_set::tree_propagate_surface_velocity(int num)
{
	if(dimension==2)
	{
		tree_propagate_velocity_fluid_points(num);

		double delta = MAXSIZE * pow(2,-max_tree_depth);

		pointphi_2d *temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->tempphi = temp->phi;
			temp->phi = temp->tempphi + temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = 3./4. * temp->tempphi + 1./4. * temp->phi + 1./4. * temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = 1./3. * temp->tempphi + 2./3. * temp->phi + 2./3. * temp->K1;
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_2d, false);

		/*pointphi_2d *temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			double xhat = temp->x - time_step/2. * external_velocity_x_2d(temp->x, temp->y, num);
			double yhat = temp->y - time_step/2. * external_velocity_y_2d(temp->x, temp->y, num);
			double xd = temp->x - time_step * external_velocity_x_2d(xhat, yhat, num);
			double yd = temp->y - time_step * external_velocity_y_2d(xhat, yhat, num);
			temp->tempphi = level_computephi(xd,yd);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = temp->tempphi;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			quadtree_coeffs_saving(tree_grid_2d, temp);
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_2d, false);*/
	}

	if(dimension==3)
	{
		tree_propagate_velocity_fluid_points(num);

		double delta = MAXSIZE * pow(2,-max_tree_depth);

		pointphi_3d *temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->tempphi = temp->phi;
			temp->phi = temp->tempphi + temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = 3./4. * temp->tempphi + 1./4. * temp->phi + 1./4. * temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_velocity(temp, delta, num);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = 1./3. * temp->tempphi + 2./3. * temp->phi + 2./3. * temp->K1;
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_3d, false);

		/*pointphi_3d *temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			double xhat = temp->x - time_step/2. * external_velocity_x_3d(temp->x, temp->y, temp->z, num);
			double yhat = temp->y - time_step/2. * external_velocity_y_3d(temp->x, temp->y, temp->z, num);
			double zhat = temp->z - time_step/2. * external_velocity_z_3d(temp->x, temp->y, temp->z, num);
			double xd = temp->x - time_step * external_velocity_x_3d(xhat, yhat, zhat, num);
			double yd = temp->y - time_step * external_velocity_y_3d(xhat, yhat, zhat, num);
			double zd = temp->z - time_step * external_velocity_z_3d(xhat, yhat, zhat, num);
			temp->tempphi = level_computephi(xd,yd,zd);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = temp->tempphi;
			temp = temp->next;
		}

		iter_reinitial_scheme();
		generate_tree(tree_grid_3d, false);*/
	}

}


double convection_normal(pointphi_2d *temp, double delta)
{
	double phi_xplus = quadtree_weno(temp,0,delta);
	double phi_xminus = quadtree_weno(temp,1,delta);
	double phi_yplus = quadtree_weno(temp,2,delta);
	double phi_yminus = quadtree_weno(temp,3,delta);

	double phi_x; double phi_y;

	if(phi_xplus>=0 && phi_xminus>=0) phi_x = phi_xminus;
	else if(phi_xplus<=0 && phi_xminus<=0) phi_x = phi_xplus;
	else if(phi_xplus>0 && phi_xminus<0) phi_x = 0;
	else if(-(abs(phi_xplus)-abs(phi_xminus)) > 0) phi_x = phi_xminus;
	else phi_x = phi_xplus;
	
	if(phi_yplus>=0 && phi_yminus>=0) phi_y = phi_yminus;
	else if(phi_yplus<=0 && phi_yminus<=0) phi_y = phi_yplus;
	else if(phi_yplus>0 && phi_yminus<0) phi_y = 0;
	else if(-(abs(phi_yplus)-abs(phi_yminus)) > 0) phi_y = phi_yminus;
	else phi_y = phi_yplus;

	return -temp->b_rate * sqrt(phi_x*phi_x+phi_y*phi_y);
}

double convection_normal(pointphi_3d *temp, double delta)
{
	double phi_xplus = octree_weno(temp,0,delta);
	double phi_xminus = octree_weno(temp,1,delta);
	double phi_yplus = octree_weno(temp,2,delta);
	double phi_yminus = octree_weno(temp,3,delta);
	double phi_zplus = octree_weno(temp,4,delta);
	double phi_zminus = octree_weno(temp,5,delta);

	double phi_x; double phi_y; double phi_z;

	if(phi_xplus>=0 && phi_xminus>=0) phi_x = phi_xminus;
	else if(phi_xplus<=0 && phi_xminus<=0) phi_x = phi_xplus;
	else if(phi_xplus>0 && phi_xminus<0) phi_x = 0;
	else if(-(abs(phi_xplus)-abs(phi_xminus)) > 0) phi_x = phi_xminus;
	else phi_x = phi_xplus;
	
	if(phi_yplus>=0 && phi_yminus>=0) phi_y = phi_yminus;
	else if(phi_yplus<=0 && phi_yminus<=0) phi_y = phi_yplus;
	else if(phi_yplus>0 && phi_yminus<0) phi_y = 0;
	else if(-(abs(phi_yplus)-abs(phi_yminus)) > 0) phi_y = phi_yminus;
	else phi_y = phi_yplus;

	if(phi_zplus>=0 && phi_zminus>=0) phi_z = phi_zminus;
	else if(phi_zplus<=0 && phi_zminus<=0) phi_z = phi_zplus;
	else if(phi_zplus>0 && phi_zminus<0) phi_z = 0;
	else if(-(abs(phi_zplus)-abs(phi_zminus)) > 0) phi_z = phi_zminus;
	else phi_z = phi_zplus;

	return -temp->b_rate * sqrt(phi_x*phi_x+phi_y*phi_y+phi_z*phi_z);
}

void level_set::tree_propagate_surface(double *burning_rate)
{
	if(dimension==2)
	{
		tree_propagate_fluid_points(burning_rate);

		tree_extrapolation(burning_rate);

		double delta = pow(2,-max_tree_depth);
		
		pointphi_2d *temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->tempphi = temp->phi;
			temp->phi = temp->tempphi + temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = 3./4. * temp->tempphi + 1./4. * temp->phi + 1./4. * temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = 1./3. * temp->tempphi + 2./3. * temp->phi + 2./3. * temp->K1;
			temp = temp->next;
		}
	
		iter_reinitial_scheme();
		generate_tree(tree_grid_2d, false);


		/*pointphi_2d *temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			//temp->tempphi = temp->phi + time_step * tree_hg(temp, false,0,0,0,0);
			temp->tempphi = temp->phi - time_step * tree_hg(temp, true,0,0,0,0);
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
			quadtree_coeffs_saving(tree_grid_2d, temp);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			//temp->tempphi = temp->phi + time_step * tree_hg(temp, false,0,0,0,0);
			temp->tempphi = temp->phi - time_step * tree_hg(temp, true,0,0,0,0);
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			temp->phi = (temp->tempphi + temp->tempphi2)/2.;
			temp = temp->next;
		}

		temp = tree_point_phi_2d_start;
		while(temp!=NULL)
		{
			quadtree_coeffs_saving(tree_grid_2d, temp);
			temp = temp->next;
		}
	
		//iter_reinitial_scheme();
		generate_tree(tree_grid_2d, false);*/
	}

	if(dimension==3)
	{
		tree_propagate_fluid_points(burning_rate);

		tree_extrapolation(burning_rate);

		double delta = pow(2,-max_tree_depth);
		
		pointphi_3d *temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->tempphi = temp->phi;
			temp->phi = temp->tempphi + temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = 3./4. * temp->tempphi + 1./4. * temp->phi + 1./4. * temp->K1;
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->K1 = time_step * convection_normal(temp, delta);
			temp = temp->next;
		}

		temp = tree_point_phi_3d_start;
		while(temp!=NULL)
		{
			temp->phi = 1./3. * temp->tempphi + 2./3. * temp->phi + 2./3. * temp->K1;
			temp = temp->next;
		}
	
		iter_reinitial_scheme();
		generate_tree(tree_grid_3d, false);
	}

}
