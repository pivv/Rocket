#include "levelsets_3propa.h"

void level_set::signphi()
{
	if(dimension ==2)
	{
		pointphi_2d *temp = tree_point_phi_2d_start;
	
		while(temp!=NULL)
		{
			temp->phi = sign(temp->phi);
			temp = temp->next;
		}
	}

	if(dimension == 3)
	{
		pointphi_3d *temp = tree_point_phi_3d_start;
	
		while(temp!=NULL)
		{
			temp->phi = sign(temp->phi);
			temp = temp->next;
		}
	}
}

void level_set::iter_reinitial_scheme()
{
	pointphi_2d *temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		temp->tempphi3 = temp->phi;
		if(temp->phi * temp->rightphi<0)
		{
			double c2 = 1./2. * minmod(quadtree_dxxphi(temp), quadtree_dxxphi(temp->right));
			double c1 = (temp->rightphi - temp->phi)/temp->rights;
			double c0 = (temp->rightphi + temp->phi)/2. - c2 * pow(temp->rights,2)/4.;

			if(abs(c2)<MINERROR) temp->si1 = temp->rights/2. + (-c0/c1);
			else if(temp->phi<0) temp->si1 = temp->rights/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->si1 = temp->rights/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
		}
		else temp->si1 = 0;
		if(temp->phi*temp->leftphi<0)
		{
			double c2 = 1./2. * minmod(quadtree_dxxphi(temp), quadtree_dxxphi(temp->left));
			double c1 = (temp->leftphi - temp->phi)/temp->lefts;
			double c0 = (temp->leftphi + temp->phi)/2. - c2 * pow(temp->lefts,2)/4.;

			if(abs(c2)<MINERROR) temp->si2 = temp->lefts/2. + (-c0/c1);
			else if(temp->phi<0) temp->si2 = temp->lefts/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->si2 = temp->lefts/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
		}
		else temp->si2 = 0;
		if(temp->phi*temp->topphi<0)
		{
			double c2 = 1./2. * minmod(quadtree_dyyphi(temp), quadtree_dyyphi(temp->top));
			double c1 = (temp->topphi - temp->phi)/temp->tops;
			double c0 = (temp->topphi + temp->phi)/2. - c2 * pow(temp->tops,2)/4.;

			if(abs(c2)<MINERROR) temp->sj1 = temp->tops/2. + (-c0/c1);
			else if(temp->phi<0) temp->sj1 = temp->tops/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->sj1 = temp->tops/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
		}
		else temp->sj1 = 0;
		if(temp->phi*temp->bottomphi<0)
		{
			double c2 = 1./2. * minmod(quadtree_dyyphi(temp), quadtree_dyyphi(temp->bottom));
			double c1 = (temp->bottomphi - temp->phi)/temp->bottoms;
			double c0 = (temp->bottomphi + temp->phi)/2. - c2 * pow(temp->bottoms,2)/4.;

			if(abs(c2)<MINERROR) temp->sj2 = temp->bottoms/2. + (-c0/c1);
			else if(temp->phi<0) temp->sj2 = temp->bottoms/2. + (-c1 + sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
			else temp->sj2 = temp->bottoms/2. + (-c1 - sqrt(pow(c1,2) - 4.*c2*c0))/(2.*c2);
		}
		else temp->sj2 = 0;
		temp = temp->next;
	}
	for(int i=0; i<10; i++) reinitial_scheme();
}


void level_set::reinitial_scheme()
{
	pointphi_2d *temp = tree_point_phi_2d_start;
	
	while(temp!=NULL)
	{
		double s1 = temp->lefts;
		double s4 = temp->rights;
		double s3 = temp->tops;
		double s2 = temp->bottoms;

		if(s1==0) s1 = MAXSIZE*10;
		if(s2==0) s2 = MAXSIZE*10;
		if(s3==0) s3 = MAXSIZE*10;
		if(s4==0) s4 = MAXSIZE*10;

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
		quadtree_coeffs_saving(tree_grid_2d, temp);
		temp = temp->next;
	}

	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		double s1 = temp->lefts;
		double s4 = temp->rights;
		double s3 = temp->tops;
		double s2 = temp->bottoms;

		if(s1==0) s1 = MAXSIZE*10;
		if(s2==0) s2 = MAXSIZE*10;
		if(s3==0) s3 = MAXSIZE*10;
		if(s4==0) s4 = MAXSIZE*10;

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
		temp->phi = (temp->tempphi + temp->tempphi2)/2.;
		if(sign(temp->phi)!=sign(temp->tempphi3))
		{
			double e = 1;
		}
		temp = temp->next;
	}
	
	temp = tree_point_phi_2d_start;
	while(temp!=NULL)
	{
		quadtree_coeffs_saving(tree_grid_2d, temp);
		temp = temp->next;
	}
}

