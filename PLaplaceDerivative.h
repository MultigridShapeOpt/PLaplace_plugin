/* Author: Jose Pinzon
   Source: https://github.com/MultigridShapeOpt
  *
  * This file is a part of the PLaplaceOptim UG4 plugin under development at 
  * the Research Group Approximation and Optimization, Hamburg University
  * and as part of the project SENSUS (LFF-GK11).
  *
  * This library is free software; you can redistribute it and/or
  * modify it under the terms of the GNU Lesser General Public
  * License as published by the Free Software Foundation; either
  * version 2.1 of the License, or (at your option) any later version.
  *
  * This library is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  * Lesser General Public License for more details.
*/
 
#ifndef P_LAPLACE_DERIVATIVE
#define P_LAPLACE_DERIVATIVE

#include <stdio.h>
#include <string>
#include "common/math/ugmath.h"
#include "bridge/bridge.h"

#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"


namespace ug{
	namespace PLaplacian{
		
		template<typename TDomain>
		class PLaplaceDerivative: public IElemDisc<TDomain>
		{
			//CONSTRUCTOR's
			public:
				PLaplaceDerivative(const char* functions, const char* subsets);
				PLaplaceDerivative(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
				
			protected:
			///	Base class type, why do we have to add it?
			typedef IElemDisc<TDomain> base_type;
				
			///Also self-type is defined
			typedef PLaplaceDerivative<TDomain> this_type;
			//PUBLIC MEMBERS:
			public:
				///	quadrature order
				bool m_bQuadOrderUserDef;
				int m_quadOrder;

				///	current shape function set
				LFEID m_lfeID;
				///	current order of disc scheme
				int m_order;
				///	current element
				GridObject* m_pElem;
				
				///	World dimension
				static const int dim = base_type::dim;
			protected:
				
				number m_p;
				number m_epsilon;
				number m_step_length;
				number m_control;

				number m_lambda_vol;
				MathVector<dim> m_vLambda_bar;
				
				number m_volume_defect;
				MathVector<dim> m_vBarycenterDefect;
				
				number m_nu_geom;
				number m_stabilization_scale;
				
				number m_elliptic_epsilon;
			

			public:
				void set_stabilization_scale(number v)
				{
					this->m_stabilization_scale=v;
				}
				void set_lambda_vol(number value){
					this-> m_lambda_vol = value;
				}
				void set_lambda_barycenter(number X, number Y, number Z){
					this-> m_vLambda_bar[0]=X;
					this-> m_vLambda_bar[1]=Y;
					
					if(this->dim ==  3){
						this-> m_vLambda_bar[2]=Z;
					}
				}
				void set_p(number value){
					this-> m_p = value;
				}
				void set_epsilon(number value){
					this-> m_epsilon = value;
				}
				void set_elliptic_epsilon(number value){
					this-> m_elliptic_epsilon = value;
				}
				void set_step_length(number value){
					this-> m_step_length = value;
				}
				void set_control(number value){
					this-> m_control = value;
				}
				
				
			public:

				void set_quad_order(size_t order);
				///	type of trial space for each function used
				virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
				///Functions for imports
				
			protected:
				///	register utils
			///	\{
				///	sets the requested assembling routines
				void set_assemble_funcs();

				void register_all_funcs(int order, int quadOrder);
				template <typename TElem, typename TFEGeom>
				void register_fe_func();
			/// \}

			//PUBLIC ASSEMBLY FUNCTIONS AND PREP FUNCTIONS, these are copy/pasted from navier_stokes_fe.h
			public:		
				
				template<typename TElem, typename TFEGeom>
				void prep_elem_loop(const ReferenceObjectID roid, const int si);

				template<typename TElem, typename TFEGeom>
				void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

				template<typename TElem, typename TFEGeom>
				void fsh_elem_loop();

				template<typename TElem, typename TFEGeom>
				void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
				
				template<typename TElem, typename TFEGeom>
				void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);
				
				template<typename TElem, typename TFEGeom>
				void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
				
				//Useless here

				template<typename TElem, typename TFEGeom>
				void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);		

				template<typename TElem, typename TFEGeom>
				void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

			protected:
				void IdentityDeformationGradient(MathMatrix<dim, dim>& IdGradU, const size_t ip);
				void DeformationVector(MathVector<dim>& W, const size_t ip);
				void DeformationGradient(MathMatrix<dim, dim>& GradW, const size_t ip);

				number heaviside(number x){
					//number diff=x-m_p;
					if(x > 0.0){
						return 1.0;
					}
					else {
						return 0.0;
					}
					return 0.0;
				}
				
				DataImport<MathVector<dim>,dim> m_imDeformationVectord1;//row1 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord2;//row2 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord3;//row3 Deformation Gradient
				
				DataImport<number, dim> m_imDeformationd1;//deformation on x-direction
				DataImport<number, dim> m_imDeformationd2;//deformation on y-direction
				DataImport<number, dim> m_imDeformationd3;//deformation on z-direction
				 
			public:
	
				// Deformation values to form deformation vector
				void set_deformation_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d1(number val);
				void set_deformation_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d2(number val);
				void set_deformation_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_deformation_d3(number val);
				// Deformation vectors for DeformationGradient and Deformation 
				void set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d1(number val);
				void set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d2(number val);
				void set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
				void set_deformation_vector_d3(number val);	
				
		};//end PLaplaceDerivative
	}//end namespace FluidOpt
}//end namespace ug
#endif /*P_LAPLACE_DERIVATIVE*/
