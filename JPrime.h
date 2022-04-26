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
 
#ifndef JPRIME_RHS
#define JPRIME_RHS

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
		class JPrime: public IElemDisc<TDomain>
		{
			//CONSTRUCTOR's
			public:
				JPrime(const char* functions, const char* subsets);
				JPrime(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
				
			protected:
			///	Base class type, why do we have to add it?
			typedef IElemDisc<TDomain> base_type;
				
			///Also self-type is defined
			typedef JPrime<TDomain> this_type;
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
				
				number m_step_length;
				number m_stabilization_scale;
			public:
				void set_stabilization_scale(number v)
				{
					this->m_stabilization_scale=v;
				}


				
				
			public:

				void set_quad_order(size_t order);
				///	type of trial space for each function used
				virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);
				///Functions for imports
				void set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user);
				void set_kinematic_viscosity(number val);
				#ifdef UG_FOR_LUA
				void set_kinematic_viscosity(const char* fctName);
				void set_kinematic_viscosity(LuaFunctionHandle fct);
				#endif
				
				//void set_nu_extension(number n);
			protected:
				///	register utils
			///	\{
				///	sets the requested assembling routines
				void set_assemble_funcs();

				void register_all_funcs(int order, int quadOrder);
				template <typename TElem, typename TFEGeom>
				void register_fe_func();
			/// \}
			///	Data import for kinematic viscosity
				DataImport<number, dim> m_imKinViscosity;

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
				void VelocityGradient(MathMatrix<dim, dim>& GradV, const size_t ip);
				void VelocityVector(MathVector<dim>& V, const size_t ip);
				void AdjointVelocityGradient(MathMatrix<dim, dim>& GradQ, const size_t ip);
				void AdjointVelocityVector(MathVector<dim>& Q, const size_t ip);

				
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd1;//row1 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd2;//row2 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imVelocityGradientd3;//row3 Velocity Gradient
				
				DataImport<number, dim> m_imVelocityd1;//flow on x-direction
				DataImport<number, dim> m_imVelocityd2;//flow on y-direction
				DataImport<number, dim> m_imVelocityd3;//flow on z-direction
				
				DataImport<number, dim> m_imAdjointVelocityd1;//flow on x-direction
				DataImport<number, dim> m_imAdjointVelocityd2;//flow on y-direction
				DataImport<number, dim> m_imAdjointVelocityd3;//flow on z-direction
				
				DataImport<number, dim> m_imPressure;//pressure value
				DataImport<number, dim> m_imAdjointPressure;//adjointpressure value
				
				
				//Adjoint Velocity vectors Q
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd1;//row1 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd2;//row2 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd3;//row3 Velocity Gradient

				
			public:
				//TODO:set default values, maybe in constructor?
				void set_pressure(SmartPtr<CplUserData<number, dim> > user);
				void set_pressure(number val);
				
				
				void set_adjoint_pressure(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_pressure(number val);
				
				//Velocity values to form vector
				void set_velocity_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d1(number val);
				void set_velocity_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d2(number val);
				void set_velocity_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_velocity_d3(number val);
				//Velocity vectors for VelocityGradient
				void set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d1(number val);
				void set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d2(number val);
				void set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_velocity_vector_d3(number val);
				// Adjoint Velocities Q gradient for AdjointVelocityGradient calculation
				void set_adjoint_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_adjoint_velocity_vector_d1(number val);
				void set_adjoint_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_adjoint_velocity_vector_d2(number val);
				void set_adjoint_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user);
				void set_adjoint_velocity_vector_d3(number val);
				//Adjoint Velocity Q values to form vector
				void set_adjoint_velocity_d1(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_velocity_d1(number val);
				void set_adjoint_velocity_d2(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_velocity_d2(number val);
				void set_adjoint_velocity_d3(SmartPtr<CplUserData<number, dim> > user);
				void set_adjoint_velocity_d3(number val);
			
				
		};//end PLaplaceDerivative
	}//end namespace FluidOpt
}//end namespace ug
#endif /*P_LAPLACE_DERIVATIVE*/
