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
 
#ifndef LARGE_PROBLEM_RHS
#define LARGE_PROBLEM_RHS

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
		class LargeProblemRHS: public IElemDisc<TDomain>
		{
			//CONSTRUCTOR's
			public:
				LargeProblemRHS(const char* functions, const char* subsets);
				LargeProblemRHS(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);
				
			protected:
			///	Base class type, why do we have to add it?
			typedef IElemDisc<TDomain> base_type;
				
			///Also self-type is defined
			typedef LargeProblemRHS<TDomain> this_type;
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
				
				number m_nu_geom;
				number m_stabilization_scale;
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
				void set_step_length(number value){
					this-> m_step_length = value;
				}
				void set_control(number value){
					this -> m_control = value;
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
				void IdentityDeformationGradient(MathMatrix<dim, dim>& IdGradU, const size_t ip);
				void DeformationVector(MathVector<dim>& W, const size_t ip);
				void DeformationGradient(MathMatrix<dim, dim>& GradW, const size_t ip);
				void VelocityGradient(MathMatrix<dim, dim>& GradV, const size_t ip);
				void VelocityVector(MathVector<dim>& V, const size_t ip);
				void AdjointVelocityGradient(MathMatrix<dim, dim>& GradQ, const size_t ip);
				void AdjointVelocityVector(MathVector<dim>& Q, const size_t ip);

				
				DataImport<MathVector<dim>,dim> m_imDeformationVectord1;//row1 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord2;//row2 Deformation Gradient
				DataImport<MathVector<dim>,dim> m_imDeformationVectord3;//row3 Deformation Gradient
				
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
				
				DataImport<number, dim> m_imDeformationd1;//deformation on x-direction
				DataImport<number, dim> m_imDeformationd2;//deformation on y-direction
				DataImport<number, dim> m_imDeformationd3;//deformation on z-direction
				
				//Adjoint Velocity vectors Q
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd1;//row1 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd2;//row2 Velocity Gradient
				DataImport<MathVector<dim>,dim> m_imAdjointVelocityGradientd3;//row3 Velocity Gradient
				
				number heaviside(number x){

					if(x > 0.0){
						return 1.0;
					}
					else {
						return 0.0;
					}
					return 0.0;
				} 
				
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
				
				void set_multiplier_vol(number val){
					this-> m_multiplier_vol=val;
				}
				void set_multiplier_bx(number val){
					this-> m_multiplier_bx=val;
				}
				void set_multiplier_by(number val){
					this-> m_multiplier_by=val;
				}
				void set_multiplier_bz(number val){
					this-> m_multiplier_bz=val;
				}
			protected:
				number m_multiplier_vol;
				number m_multiplier_bx;
				number m_multiplier_by;
				number m_multiplier_bz;
			
				
		};//end PLaplaceDerivative
	}//end namespace FluidOpt
}//end namespace ug
#endif /*P_LAPLACE_DERIVATIVE*/
