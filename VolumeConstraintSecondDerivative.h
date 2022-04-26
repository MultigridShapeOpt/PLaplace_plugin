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
 
#ifndef VOLUME_SECOND_DERIVATIVE
#define VOLUME_SECOND_DERIVATIVE

#include <stdio.h>
#include <string>

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
		class VolumeConstraintSecondDerivative: public IElemDisc<TDomain>
		{
			//CONSTRUCTOR's
			public:
				VolumeConstraintSecondDerivative(const char* functions, const char* subsets);
				VolumeConstraintSecondDerivative(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
				
			protected:
			///	Base class type, why do we have to add it?
			typedef IElemDisc<TDomain> base_type;
				
			///Also self-type is defined
			typedef VolumeConstraintSecondDerivative<TDomain> this_type;
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

				
			public:

				void set_quad_order(size_t order);
				///	type of trial space for each function used
				virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

			protected:
				///	register utils
			///	\{
				///	sets the requested assembling routines
				void set_assemble_funcs();

				void register_all_funcs(int order, int quadOrder);
				template <typename TElem, typename TFEGeom>
				void register_fe_func();
				
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
				
				void set_multiplier(number val){
					this->m_multiplier=val;
				}
			
				number m_multiplier;
				
				
		};//end VolumeConstraintSecondDerivative
	}//end namespace PLaplacian
}//end namespace ug
#endif /*P_LAPLACE_DERIVATIVE*/
