#ifndef ELEMENT_H
#define ELEMENT_H

class Element {
public:
	int id;   /**< Element identification. */
	int q;    /**< Quadrature order. */
	int L;    /**< Polynomial order. */
	int *Z;   /**< Global numbering. */

	double x[2]; /**< Left \f$ x_{e}^{L} \f$ and right \f$ x_{e}^{R} \f$ axial coordinates of each element. */
	double rx; /**< Inverse of the Jacobian "jac". */
	double jac; /**< Jacobian of the elemenetal mapping onto the standard element.
				\f$J_{e} = \frac{1}{2} (x_{e}^{R} - x_{e}^{L})\f$*/
	double *beta; /**< material properties "Beta" in each quadrature point. */
	double *Ao; /**< initial area "Ao" in each quadrature point. */
	double *Hdo;

	double *h; /**< Solution in the quadrature points of the physical space. */
	double *hj; /**< Solution in the projected space. */

	double *Xgeom;
	double *Ygeom;

	Element *next; /**< Pointer to the next element of the domain. */

	double Get_P  (int i);
	void   Get_P  (double *p);
	double Get_P  (int i, double A);
	double PtoA   (int i, double p);
	double Get_W  (int i);
	void   Get_W  (double *p);
	double Get_W  (int i, double A);
	double Get_c  (int i);
	double Get_c  (int i, double A);
	void   Get_c  (double *p);
	double Get_co (int i);
	void   Get_co (double *p);

private:
	double eval_A(double p, double ao, double beta);
	double eval_P(double A, double ao, double beta);
	double eval_c(double A, double beta);
	double eval_W(double A, double ao, double beta);
	/*double eval_BC_A(double A_star, double A);
	double eval_A_from_P(double p, double po, double beta, double Ao);*/
};

#endif
