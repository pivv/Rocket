#include "trees.h"

void quadtree_coeffs_saving(struct quadtree *start, struct pointphi_2d *point)		// phi0, phi1, s1, phi4, s4, phi3, s3, phi2, s2
{
	if(point->left!=NULL)
	{
		point->leftphi = point->left->phi;
		point->lefts = point->x - point->left->x;
	}
	else if(point->x==0)
	{
		point->leftphi = 0;
		point->lefts = 0;
	}
	
	if(point->right!=NULL)
	{
		point->rightphi = point->right->phi;
		point->rights = point->right->x - point->x;
	}
	else if(point->x==MAXSIZE)
	{
		point->rightphi = 0;
		point->rights = 0;
	}
	
	if(point->top!=NULL)
	{
		point->topphi = point->top->phi;
		point->tops = point->top->y - point->y;
	}
	else if(point->y==MAXSIZE)
	{
		point->topphi = 0;
		point->tops = 0;
	}
	
	if(point->bottom!=NULL)
	{
		point->bottomphi = point->bottom->phi;
		point->bottoms = point->y - point->bottom->y;
	}
	else if(point->y==0)
	{
		point->bottomphi = 0;
		point->bottoms = 0;
	}

	if(point->left==NULL && point->x!=0)
	{
		struct quadtree *left = find_left_top_face(start, point->x, point->y);
		double phi5 = left->phi_lefttop->phi;
		double phi6 = left->phi_leftbottom->phi;
		double s5 = left->phi_lefttop->y - point->y;
		double s6 = point->y - left->phi_leftbottom->y;
		
		point->leftphi = (phi5*s6 + phi6*s5)/(s5+s6) -
			s5*s6/(point->bottoms+point->tops) * ((point->bottomphi-point->phi)/point->bottoms + (point->topphi-point->phi)/point->tops);
		point->lefts = point->x - left->phi_lefttop->x;
	}
	else if(point->right==NULL && point->x!=MAXSIZE)
	{
		struct quadtree *right = find_right_bottom_face(start, point->x, point->y);
		double phi5 = right->phi_righttop->phi;
		double phi6 = right->phi_rightbottom->phi;
		double s5 = right->phi_righttop->y - point->y;
		double s6 = point->y - right->phi_rightbottom->y;

		point->rightphi = (phi5*s6 + phi6*s5)/(s5+s6) -
			s5*s6/(point->bottoms+point->tops) * ((point->bottomphi-point->phi)/point->bottoms + (point->topphi-point->phi)/point->tops);
		point->rights = right->phi_righttop->x - point->x;
	}
	else if(point->top==NULL && point->y!=MAXSIZE)
	{
		struct quadtree *top = find_left_top_face(start, point->x, point->y);
		double phi5 = top->phi_lefttop->phi;
		double phi6 = top->phi_righttop->phi;
		double s5 = point->x - top->phi_lefttop->x;
		double s6 = top->phi_righttop->x - point->x;

		point->topphi = (phi5*s6 + phi6*s5)/(s5+s6) -
			s5*s6/(point->lefts+point->rights) * ((point->leftphi-point->phi)/point->lefts + (point->rightphi-point->phi)/point->rights);
		point->tops = top->phi_lefttop->y - point->y;
	}
	else if(point->bottom==NULL && point->y!=0)
	{
		struct quadtree *bottom = find_right_bottom_face(start, point->x, point->y);
		double phi5 = bottom->phi_leftbottom->phi;
		double phi6 = bottom->phi_rightbottom->phi;
		double s5 = point->x - bottom->phi_leftbottom->x;
		double s6 = bottom->phi_rightbottom->x - point->x;

		point->bottomphi = (phi5*s6 + phi6*s5)/(s5+s6) -
			s5*s6/(point->lefts+point->rights) * ((point->leftphi-point->phi)/point->lefts + (point->rightphi-point->phi)/point->rights);
		point->bottoms = point->y - bottom->phi_leftbottom->y;
	}
}

inline double quadtree_dxphi(struct pointphi_2d *point)
{
	if(point==NULL) return 0;
	if(point->lefts==0 || point->rights==0) return 0;
	return (point->rightphi - point->phi)/point->rights * point->lefts/(point->lefts + point->rights) +
		(point->phi - point->leftphi)/point->lefts * point->rights/(point->lefts + point->rights);
}

inline double quadtree_dxxphi(struct pointphi_2d *point)
{
	if(point==NULL) return 0;
	if(point->lefts==0 || point->rights==0) return 0;
	return (point->rightphi - point->phi)/point->rights * 2./(point->lefts + point->rights) -
		(point->phi - point->leftphi)/point->lefts * 2./(point->lefts + point->rights);
}

inline double quadtree_dxplusphi(struct pointphi_2d *point, double si)
{
	if(point==NULL) return 0;
	if(point->rights==0) return 0;
	if(si==0) return (point->rightphi - point->phi)/point->rights
		- point->rights/2. * minmod(quadtree_dxxphi(point), quadtree_dxxphi(point->right));
	else return (0 - point->phi)/si - si/2. * minmod(quadtree_dxxphi(point), quadtree_dxxphi(point->right));
}

inline double quadtree_dxminusphi(struct pointphi_2d *point, double si)
{
	if(point==NULL) return 0;
	if(point->lefts==0) return 0;
	if(si==0) return (point->phi - point->leftphi)/point->lefts
		+ point->lefts/2. * minmod(quadtree_dxxphi(point), quadtree_dxxphi(point->left));
	else return (point->phi - 0)/si + si/2. * minmod(quadtree_dxxphi(point), quadtree_dxxphi(point->left));
}

inline double quadtree_dyphi(struct pointphi_2d *point)
{
	if(point==NULL) return 0;
	if(point->tops==0 || point->bottoms==0) return 0;
	return (point->topphi - point->phi)/point->tops * point->bottoms/(point->tops + point->bottoms) +
		(point->phi - point->bottomphi)/point->bottoms * point->tops/(point->tops + point->bottoms);
}


inline double quadtree_dyyphi(struct pointphi_2d *point)
{
	if(point==NULL) return 0;
	if(point->tops==0 || point->bottoms==0) return 0;
	return (point->topphi - point->phi)/point->tops * 2./(point->tops + point->bottoms) -
		(point->phi - point->bottomphi)/point->bottoms * 2./(point->tops + point->bottoms);
}

inline double quadtree_dyplusphi(struct pointphi_2d *point, double si)
{
	if(point==NULL) return 0;
	if(point->tops==0) return 0;
	if(si==0) return (point->topphi - point->phi)/point->tops
		- point->tops/2. * minmod(quadtree_dyyphi(point), quadtree_dyyphi(point->top));
	else return (0 - point->phi)/si - si/2. * minmod(quadtree_dyyphi(point), quadtree_dyyphi(point->top));
}

inline double quadtree_dyminusphi(struct pointphi_2d *point, double si)
{
	if(point==NULL) return 0;
	if(point->bottoms==0) return 0;
	if(si==0) return (point->phi - point->bottomphi)/point->bottoms
		- point->bottoms/2. * minmod(quadtree_dyyphi(point), quadtree_dyyphi(point->bottom));
	else return (point->phi - 0)/si - si/2. * minmod(quadtree_dyyphi(point), quadtree_dyyphi(point->bottom));
}




inline double octree_dxphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->left==NULL || point->right==NULL) return 0;

	double lefts = point->x - point->left->x;
	double rights = point->right->x - point->x;
	return (point->right->phi - point->phi)/rights * lefts/(lefts + rights) +
		(point->phi - point->left->phi)/lefts * rights/(lefts + rights);
}

inline double octree_dxxphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->left==NULL || point->right==NULL) return 0;

	double lefts = point->x - point->left->x;
	double rights = point->right->x - point->x;
	return (point->right->phi - point->phi)/rights * 2./(lefts + rights) -
		(point->phi - point->left->phi)/lefts * 2./(lefts + rights);
}

inline double octree_dyphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->bottom==NULL || point->top==NULL) return 0;

	double bottoms = point->x - point->bottom->x;
	double tops = point->top->x - point->x;
	return (point->top->phi - point->phi)/tops * bottoms/(bottoms + tops) -
		(point->phi - point->bottom->phi)/bottoms * tops/(bottoms + tops);
}

inline double octree_dyyphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->bottom==NULL || point->top==NULL) return 0;

	double bottoms = point->y - point->bottom->y;
	double tops = point->top->y - point->y;
	return (point->top->phi - point->phi)/tops * 2./(bottoms + tops) -
		(point->phi - point->bottom->phi)/bottoms * 2./(bottoms + tops);
}

inline double octree_dzphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->back==NULL || point->front==NULL) return 0;

	double backs = point->z - point->back->z;
	double fronts = point->front->z - point->z;
	return (point->front->phi - point->phi)/fronts * backs/(backs + fronts) -
		(point->phi - point->back->phi)/backs * fronts/(backs + fronts);
}

inline double octree_dzzphi(struct pointphi_3d *point)
{
	if(point==NULL) return 0;
	if(point->back==NULL || point->front==NULL) return 0;

	double backs = point->z - point->back->z;
	double fronts = point->front->z - point->z;
	return (point->front->phi - point->phi)/fronts * 2./(backs + fronts) -
		(point->phi - point->back->phi)/backs * 2./(backs + fronts);
}

double tree_hg(struct pointphi_2d *point, bool pm, double si1, double si2, double sj1, double sj2)
{
	double a;
	double b;
	double c;
	double d;

	if(pm==true)
	{
		a = quadtree_dxplusphi(point,si1);
		b = quadtree_dxminusphi(point,si2);
		c = quadtree_dyplusphi(point,sj1);
		d = quadtree_dyminusphi(point,sj2);
	}
	else
	{
		a = quadtree_dxminusphi(point,si2);
		b = quadtree_dxplusphi(point,si1);
		c = quadtree_dyminusphi(point,sj2);
		d = quadtree_dyplusphi(point,sj1);
	}

	return sqrt(max(pow(abs(min(a,0)),2),pow(abs(max(b,0)),2)) + max(pow(abs(min(c,0)),2),pow(abs(max(d,0)),2)));
}




double computephi(struct quadtree *current, double x, double y)
{
	double x1 = current->phi_lefttop->x;
	double x2 = current->phi_righttop->x;
	double y1 = current->phi_leftbottom->y;
	double y2 = current->phi_lefttop->y;
	double l = current->length;

	double a = abs(quadtree_dxxphi(current->phi_lefttop));
	double b = abs(quadtree_dxxphi(current->phi_leftbottom));
	double c = abs(quadtree_dxxphi(current->phi_righttop));
	double d = abs(quadtree_dxxphi(current->phi_rightbottom));

	double phixx = min(min(min(a, b), c), d);

	a = abs(quadtree_dyyphi(current->phi_lefttop));
	b = abs(quadtree_dyyphi(current->phi_leftbottom));
	c = abs(quadtree_dyyphi(current->phi_righttop));
	d = abs(quadtree_dyyphi(current->phi_rightbottom));

	double phiyy = min(min(min(a, b), c), d);

	return current->phi_leftbottom->phi*(x2-x)/l*(y2-y)/l + current->phi_lefttop->phi*(x2-x)/l*(y-y1)/l
		+ current->phi_rightbottom->phi*(x-x1)/l*(y2-y)/l + current->phi_righttop->phi*(x-x1)/l*(y-y1)/l
		- phixx*(x-x1)*(x2-x)/2. - phiyy*(y-y1)*(y2-y)/2.;
}

double computephi(struct octree *current, double x, double y, double z)
{
	double x1 = current->phi_leftbottomback->x;
	double x2 = current->phi_rightbottomback->x;
	double y1 = current->phi_leftbottomback->y;
	double y2 = current->phi_lefttopback->y;
	double z1 = current->phi_leftbottomback->z;
	double z2 = current->phi_leftbottomfront->z;
	double l = current->length;
	
	/*double phixx = abs(octree_dxxphi(current->phi_leftbottomback));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_lefttopback)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_rightbottomback)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_righttopback)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_leftbottomfront)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_lefttopfront)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_rightbottomfront)));
	phixx = min(phixx, abs(octree_dxxphi(current->phi_righttopfront)));
	
	double phiyy = abs(octree_dyyphi(current->phi_leftbottomback));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_lefttopback)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_rightbottomback)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_righttopback)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_leftbottomfront)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_lefttopfront)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_rightbottomfront)));
	phiyy = min(phiyy, abs(octree_dyyphi(current->phi_righttopfront)));*/

	return current->phi_leftbottomback->phi*(x2-x)/l*(y2-y)/l*(z2-z)/l + current->phi_lefttopback->phi*(x2-x)/l*(y-y1)/l*(z2-z)/l
		+ current->phi_rightbottomback->phi*(x-x1)/l*(y2-y)/l*(z2-z)/l + current->phi_righttopback->phi*(x-x1)/l*(y-y1)/l*(z2-z)/l
		+ current->phi_leftbottomfront->phi*(x2-x)/l*(y2-y)/l*(z-z1)/l + current->phi_lefttopfront->phi*(x2-x)/l*(y-y1)/l*(z-z1)/l
		+ current->phi_rightbottomfront->phi*(x-x1)/l*(y2-y)/l*(z-z1)/l + current->phi_righttopfront->phi*(x-x1)/l*(y-y1)/l*(z-z1)/l;
}




/*double *quadtree_coeffs(struct quadtree *start, struct pointphi_2d *point)		// phi0, phi1, s1, phi4, s4, phi3, s3, phi2, s2
{
	if(point==NULL) return NULL;

	double *r = new double[9];
	r[0] = point->phi;
	if(point->left!=NULL)
	{
		r[1] = point->left->phi;
		r[2] = point->x - point->left->x;
	}
	else if(point->x==0)
	{
		r[1] = 0;
		r[2] = 0;
	}
	
	if(point->right!=NULL)
	{
		r[3] = point->right->phi;
		r[4] = point->right->x - point->x;
	}
	else if(point->x==MAXSIZE)
	{
		r[3] = 0;
		r[4] = 0;
	}
	
	if(point->top!=NULL)
	{
		r[5] = point->top->phi;
		r[6] = point->top->y - point->y;
	}
	else if(point->y==MAXSIZE)
	{
		r[5] = 0;
		r[6] = 0;
	}
	
	if(point->bottom!=NULL)
	{
		r[7] = point->bottom->phi;
		r[8] = point->y - point->bottom->y;
	}
	else if(point->y==0)
	{
		r[7] = 0;
		r[8] = 0;
	}

	if(point->left==NULL && point->x!=0)
	{
		struct quadtree *left = find_left_top_face(start, point->x, point->y);
		double phi5 = left->phi_lefttop->phi;
		double phi6 = left->phi_leftbottom->phi;
		double s5 = left->phi_lefttop->y - point->y;
		double s6 = point->y - left->phi_leftbottom->y;
		
		r[1] = (phi5*s6 + phi6*s5)/(s5+s6) - s5*s6/(r[8]+r[6]) * ((r[7]-r[0])/r[8] + (r[5]-r[0])/r[6]);
		r[2] = point->x - left->phi_lefttop->x;
	}
	else if(point->right==NULL && point->x!=MAXSIZE)
	{
		struct quadtree *right = find_right_bottom_face(start, point->x, point->y);
		double phi5 = right->phi_righttop->phi;
		double phi6 = right->phi_rightbottom->phi;
		double s5 = right->phi_righttop->y - point->y;
		double s6 = point->y - right->phi_rightbottom->y;

		r[3] = (phi5*s6 + phi6*s5)/(s5+s6) - s5*s6/(r[8]+r[6]) * ((r[7]-r[0])/r[8] + (r[5]-r[0])/r[6]);
		r[4] = right->phi_righttop->x - point->x;
	}
	else if(point->top==NULL && point->y!=MAXSIZE)
	{
		struct quadtree *top = find_left_top_face(start, point->x, point->y);
		double phi5 = top->phi_lefttop->phi;
		double phi6 = top->phi_righttop->phi;
		double s5 = point->x - top->phi_lefttop->x;
		double s6 = top->phi_righttop->x - point->x;

		r[5] = (phi5*s6 + phi6*s5)/(s5+s6) - s5*s6/(r[2]+r[4]) * ((r[1]-r[0])/r[2] + (r[3]-r[0])/r[4]);
		r[6] = top->phi_lefttop->y - point->y;
	}
	else if(point->bottom==NULL && point->y!=0)
	{
		struct quadtree *bottom = find_right_bottom_face(start, point->x, point->y);
		double phi5 = bottom->phi_leftbottom->phi;
		double phi6 = bottom->phi_rightbottom->phi;
		double s5 = point->x - bottom->phi_leftbottom->x;
		double s6 = bottom->phi_rightbottom->x - point->x;

		r[7] = (phi5*s6 + phi6*s5)/(s5+s6) - s5*s6/(r[2]+r[4]) * ((r[1]-r[0])/r[2] + (r[3]-r[0])/r[4]);
		r[8] = point->y - bottom->phi_leftbottom->y;
	}

	return r;
}

inline double quadtree_dxphi(double *r)
{
	if(r==NULL) return 0;
	if(r[2]==0 || r[4]==0) return 0;
	return (r[3] - r[0])/r[4] * r[2]/(r[2] + r[4]) + (r[0] - r[1])/r[2] * r[4]/(r[2] + r[4]);
}

inline double quadtree_dxxphi(double *r)
{
	if(r==NULL) return 0;
	if(r[2]==0 || r[4]==0) return 0;
	return (r[3] - r[0])/r[4] * 2./(r[2] + r[4]) - (r[0] - r[1])/r[2] * 2./(r[2] + r[4]);
}

inline double quadtree_dxplusphi(double *r, double *r2, double si)
{
	if(r==NULL) return 0;
	if(r[4]==0) return 0;
	if(si==0) return (r[3] - r[0])/r[4]; //- r[4]/2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
	else return (0 - r[0])/si - si/2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
}

inline double quadtree_dxminusphi(double *r, double *r2, double si)
{
	if(r==NULL) return 0;
	if(r[2]==0) return 0;
	if(si==0) return (r[0] - r[1])/r[2]; //+ r[2]/2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
	else return (r[0] - 0)/si + si/2. * minmod(quadtree_dxxphi(r), quadtree_dxxphi(r2));
}

inline double quadtree_dyphi(double *r)
{
	if(r==NULL) return 0;
	if(r[6]==0 || r[8]==0) return 0;
	return (r[5] - r[0])/r[6] * r[8]/(r[6] + r[8]) + (r[0] - r[7])/r[8] * r[6]/(r[6] + r[8]);
}


inline double quadtree_dyyphi(double *r)
{
	if(r==NULL) return 0;
	if(r[6]==0 || r[8]==0) return 0;
	return (r[5] - r[0])/r[6] * 2./(r[6] + r[8]) - (r[0] - r[7])/r[8] * 2./(r[6] + r[8]);
}

inline double quadtree_dyplusphi(double *r, double *r2, double si)
{
	if(r==NULL) return 0;
	if(r[6]==0) return 0;
	if(si==0) return (r[5] - r[0])/r[6]; //- r[6]/2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
	else return (0 - r[0])/si - si/2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
}

inline double quadtree_dyminusphi(double *r, double *r2, double si)
{
	if(r==NULL) return 0;
	if(r[8]==0) return 0;
	if(si==0) return (r[0] - r[7])/r[8]; //- r[8]/2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
	else return (r[0] - 0)/si - si/2. * minmod(quadtree_dyyphi(r), quadtree_dyyphi(r2));
}

double computephi_old(struct quadtree *start, struct quadtree *current, double x, double y)
{
	double x1 = current->phi_lefttop->x;
	double x2 = current->phi_righttop->x;
	double y1 = current->phi_leftbottom->y;
	double y2 = current->phi_lefttop->y;
	double l = current->length;

	double *r1 = quadtree_coeffs(start, current->phi_lefttop);
	double *r2 = quadtree_coeffs(start, current->phi_leftbottom);
	double *r3 = quadtree_coeffs(start, current->phi_righttop);
	double *r4 = quadtree_coeffs(start, current->phi_rightbottom);

	double a;
	double b;
	double c;
	double d;

	if(r1==NULL) a = MAXSIZE*10;
	else a = abs(quadtree_dxxphi(r1));
	if(r2==NULL) b = MAXSIZE*10;
	else b = abs(quadtree_dxxphi(r2));
	if(r3==NULL) c = MAXSIZE*10;
	else c = abs(quadtree_dxxphi(r3));
	if(r4==NULL) d = MAXSIZE*10;
	else d = abs(quadtree_dxxphi(r4));

	double phixx = min(min(min(a, b), c), d);

	if(r1!=NULL) a = abs(quadtree_dxxphi(r1));
	if(r2!=NULL) b = abs(quadtree_dxxphi(r2));
	if(r3!=NULL) c = abs(quadtree_dxxphi(r3));
	if(r4!=NULL) d = abs(quadtree_dxxphi(r4));

	double phiyy = min(min(min(a, b), c), d);

	if(r1!=NULL) delete r1;
	if(r2!=NULL) delete r2;
	if(r3!=NULL) delete r3;
	if(r4!=NULL) delete r4;

	return current->phi_leftbottom->phi*(x2-x)/l*(y2-y)/l + current->phi_lefttop->phi*(x2-x)/l*(y-y1)/l
		+ current->phi_rightbottom->phi*(x-x1)/l*(y2-y)/l + current->phi_righttop->phi*(x-x1)/l*(y-y1)/l
		- phixx*(x-x1)/l*(x2-x)/l/2. - phiyy*(y-y1)/l*(y2-y)/l/2.;
}*/