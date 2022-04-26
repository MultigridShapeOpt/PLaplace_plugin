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

#include "JPrime.h"

// for various user data
#include "bindings/lua/lua_user_data.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"

using namespace std;
namespace ug{
namespace PLaplacian{
////////////////////////////////////////////////////////////////////////////////
		//NON-REGISTERED FUNCTIONS: THEY USE IMPORTS TO CALCULATE, NOT PUBLIC
		////////////////////////////////////////////////////////////////////////////////
		//Calculate a matrix using the imports for velocity gradient
		template<typename TDomain>
		void JPrime<TDomain>::
		VelocityGradient(MathMatrix<dim, dim>& GradV, const size_t ip)
		{
			MathVector<dim> velocity_d1=m_imVelocityGradientd1[ip];
			MathVector<dim> velocity_d2=m_imVelocityGradientd2[ip];
			
			GradV.assign(velocity_d1,0);
			GradV.assign(velocity_d2,1);
			if(this->dim == 3){
				MathVector<dim> velocity_d3=m_imVelocityGradientd3[ip];
				GradV.assign(velocity_d3,2);
			}
		}
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void JPrime<TDomain>::
		VelocityVector(MathVector<dim>& V, const size_t ip)
		{
			number v1= m_imVelocityd1[ip];
			number v2= m_imVelocityd2[ip];
			
			V[0]=v1;
			V[1]=v2;
			if(this->dim == 3){
				number v3= m_imVelocityd3[ip];
				V[2]=v3;
			}
		}

		//Calculate a matrix using the imports for adjoint velocity gradient
		template<typename TDomain>
		void JPrime<TDomain>::
		AdjointVelocityGradient(MathMatrix<dim, dim>& AdjGradV, const size_t ip)
		{
			MathVector<dim> adjoint_velocity_d1=m_imAdjointVelocityGradientd1[ip];
			MathVector<dim> adjoint_velocity_d2=m_imAdjointVelocityGradientd2[ip];
			
			AdjGradV.assign(adjoint_velocity_d1,0);
			AdjGradV.assign(adjoint_velocity_d2,1);
			if(this->dim == 3){
				MathVector<dim> adjoint_velocity_d3=m_imAdjointVelocityGradientd3[ip];
				AdjGradV.assign(adjoint_velocity_d3,2);
			}
		}
		//Calculates the adjoint velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void JPrime<TDomain>::
		AdjointVelocityVector(MathVector<dim>& Q, const size_t ip)
		{
			number q1= m_imAdjointVelocityd1[ip];
			number q2= m_imAdjointVelocityd2[ip];
			
			Q[0]=q1;
			Q[1]=q2;
			if(this->dim == 3){
				number q3= m_imAdjointVelocityd3[ip];
				Q[2]=q3;
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		//	IMPORT SETTING FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////
		
		//Let's not use the LUA handlers for these imports, since we pass a CplUserData derived class directly :)
		//***************IMPORT: PRESSURE SCALAR VALUE
		template<typename TDomain>
		void JPrime<TDomain>::
		set_pressure(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imPressure.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_pressure(number val)
		{
			set_pressure(make_sp(new ConstUserNumber<dim>(val)));
		}//***************IMPORT: ADJOINT PRESSURE SCALAR VALUE
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_pressure(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointPressure.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_pressure(number val)
		{
			set_adjoint_pressure(make_sp(new ConstUserNumber<dim>(val)));
		}

		//***************IMPORT: ADJOINT VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d1(number val)
		{
			set_adjoint_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d2(number val)
		{
			set_adjoint_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_d3(number val)
		{
			set_adjoint_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: ADJOINT VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d1(number val)
		{
			set_adjoint_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d2(number val)
		{
			set_adjoint_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_adjoint_velocity_vector_d3(number val)
		{
			set_adjoint_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		//***************IMPORT: VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d1(number val)
		{
			set_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d2(number val)
		{
			set_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_d3(number val)
		{
			set_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d1(number val)
		{
			set_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d2(number val)
		{
			set_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_velocity_vector_d3(number val)
		{
			set_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		
		template<typename TDomain>
		void JPrime<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void JPrime<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void JPrime<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void JPrime<TDomain>::
		set_kinematic_viscosity(LuaFunctionHandle fct)
		{
			set_kinematic_viscosity(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif

		

		////////////////////////////////////////////////////////////////////////////////
		//	General
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		void JPrime<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void JPrime<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("JPrime: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("JPrime: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

			//	update assemble functions
			set_assemble_funcs();
		}
		template <typename TDomain>
		void
		JPrime<TDomain>::
		set_assemble_funcs()
		{
			//	set default quadrature order if not set by user
			if (!m_bQuadOrderUserDef) {
				m_quadOrder = 2 * m_order + 1;
			}
			//	set all non-set orders
			else
			{
				if (m_quadOrder < 0){
					m_quadOrder = 2 * m_order + 1;
				}
			}

			register_all_funcs(m_order, m_quadOrder);

		}
		////////////////////////////////////////////////////////////////////////////////
		// Assembling functions
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{

			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			//	prepare geometry for type and order
			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}UG_CATCH_THROW("JPrime::prep_elem_loop:"
							" Cannot update Finite Element Geometry.");

			m_imKinViscosity.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imVelocityGradientd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imVelocityGradientd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imVelocityGradientd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imVelocityd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imVelocityd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imVelocityd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imAdjointVelocityd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imAdjointVelocityd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imAdjointVelocityd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imAdjointPressure.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imPressure.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imAdjointVelocityGradientd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imAdjointVelocityGradientd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imAdjointVelocityGradientd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
		}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::fsh_elem_loop()
		{}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("JPrime::prep_elem:"
							" Cannot update Finite Element Geometry.");
			m_imKinViscosity.set_global_ips(geo.global_ips(), geo.num_ip());

			
			m_imVelocityGradientd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imVelocityGradientd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imVelocityGradientd3.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imVelocityd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imVelocityd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imVelocityd3.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imAdjointVelocityd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imAdjointVelocityd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imAdjointVelocityd3.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imAdjointPressure.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imPressure.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imAdjointVelocityGradientd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imAdjointVelocityGradientd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imAdjointVelocityGradientd3.set_global_ips(geo.global_ips(), geo.num_ip());
			
		}
		//  assemble stiffness jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	
			//	request geometry
		}
		//  assemble right-hand-side d(i,sh)
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		 
			const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);


			for (size_t ip = 0; ip < geo.num_ip(); ++ip){
								
				
				//Vectors
				MathVector<dim> vVelocity;VecSet(vVelocity,0.0);
				MathVector<dim> vAdjVelocity;VecSet(vAdjVelocity,0.0);
				
				AdjointVelocityVector(vAdjVelocity,ip);
				VelocityVector(vVelocity,ip);
				
				//Declaration of import data structures

				MathMatrix<dim, dim> AdjVelGrad;MatSet(AdjVelGrad,0.0);
				MathMatrix<dim, dim> VelGrad;MatSet(VelGrad,0.0);
				const number pressure=m_imPressure[ip];
				const number adjoint_pressure=m_imAdjointPressure[ip];
				const number viscosity = m_imKinViscosity[ip];				
				//Calculation
			
				VelocityGradient(VelGrad,ip);
				AdjointVelocityGradient(AdjVelGrad,ip);
				
				
				number TrAdjVel=Trace(AdjVelGrad);		
				number TrVel=Trace(VelGrad);

			
				for (int d1 = 0; d1 < dim; ++d1){
					for (size_t sh1 = 0; sh1 < geo.num_sh(); ++sh1){
						
						//*******************************************
						//******J'(v), as part of RHS [SENSITIVITY]**
						//*******************************************
					
						number TrDdelta_u=0.0;
						number TrAdjVelDdelta_u = 0.0;
						number TrVelDdelta_u=0.0;
						number TrDefInvDdelta_u=0.0;
						MathMatrix<dim, dim> Ddelta_u;MatSet(Ddelta_u,0.0);
						Ddelta_u.assign(geo.global_grad(ip,sh1),d1);
						TrDdelta_u=Trace(Ddelta_u);
						
							
						//
						MathMatrix<dim, dim> AdjointVelDdelta_u;MatSet(AdjointVelDdelta_u,0.0);
						MatMultiply(AdjointVelDdelta_u,AdjVelGrad,Ddelta_u);
						TrAdjVelDdelta_u = Trace(AdjointVelDdelta_u);
						
						//
						MathMatrix<dim, dim> VelocityDdelta_u;MatSet(VelocityDdelta_u,0.0);
						MatMultiply(VelocityDdelta_u,VelGrad,Ddelta_u);
						TrVelDdelta_u = Trace(VelocityDdelta_u);
												
						//P-Terms
						d(d1, sh1) -= m_step_length*TrAdjVel*TrDdelta_u*geo.weight(ip)*pressure;//p-term 1
						d(d1, sh1) += m_step_length*TrAdjVelDdelta_u*geo.weight(ip)*pressure;//p-term 2		
						//H-Terms
						d(d1, sh1) -= m_step_length*TrVel*TrDdelta_u*geo.weight(ip)*adjoint_pressure;//h-term 1
						d(d1, sh1) += m_step_length*TrVelDdelta_u*geo.weight(ip)*adjoint_pressure;//h-term 1
						
						//For Viscosity-Terms
						number contraction_velocity=MatContraction(VelocityDdelta_u,VelGrad);
						number contraction_adjoint_velocity=MatContraction(VelocityDdelta_u,AdjVelGrad);
						//For term 3

						number contraction_velocity2=MatContraction(AdjointVelDdelta_u,VelGrad);
						
						//For term 4, now 2
						MathMatrix<dim, dim> ScaledVelGrad;MatSet(ScaledVelGrad,0.0);
						MatMultiply(ScaledVelGrad,VelGrad,TrDdelta_u);
						//number contraction_scaled_velocity=MatContraction(ScaledVelGrad,VelGrad);//TODO:why was I so stupid with this term?
						number contraction_scaled_velocity = TrDdelta_u*MatContraction(VelGrad,VelGrad);
						//For term 5
						number contraction_velocity_adjointvel=MatContraction(VelGrad,AdjVelGrad);
						
						d(d1, sh1) += m_step_length*viscosity*contraction_velocity*geo.weight(ip);//Viscosity-Term 1
						d(d1, sh1) -= m_step_length*0.5*viscosity*contraction_scaled_velocity*geo.weight(ip);//Viscosity-Term 2
						d(d1, sh1) -= m_step_length*viscosity*contraction_adjoint_velocity*geo.weight(ip);//Viscosity-Term 3
						d(d1, sh1) -= m_step_length*viscosity*contraction_velocity2*geo.weight(ip);//Viscosity-Term 4
						d(d1, sh1) += m_step_length*viscosity*contraction_velocity_adjointvel*TrDdelta_u*geo.weight(ip);//Viscosity-Term 5
												
						//Terms that include the dot product with Q(Adjoint velocity)
						MathVector<dim> vTemp1;VecSet(vTemp1,0.0);
						MatVecMult(vTemp1,VelocityDdelta_u,vVelocity);
						number dot_product=VecDot(vTemp1,vAdjVelocity);//Grad_V*Grad_Test_W*V.Q
						
						MathVector<dim> vTemp2;VecSet(vTemp2,0.0);
						MatVecMult(vTemp2, VelGrad, vVelocity);
						number dot_product2=VecDot(vTemp2,vAdjVelocity);//Grad_V*V.Q
						
						//Both as a single term now
						//d(d1, sh1) -= m_step_length*(-dot_product + TrDdelta_u*dot_product2) * geo.weight(ip);//Q.terms 6-7
						//Separate as 2 terms and proper variable naming
						d(d1, sh1) -= m_step_length*  dot_product * geo.weight(ip);
						d(d1, sh1) += m_step_length* TrDdelta_u * dot_product2 * geo.weight(ip);						
													
					}
				}
				
			}// end for ip
			
	
			
		}
		//  assemble stiffness defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}//end add defect
		
		//	assemble mass jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//  assemble mass-defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}

		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		JPrime<TDomain>::
		JPrime(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			m_step_length = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'JPrime'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports
			this->register_import(m_imKinViscosity);			
			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->register_import(m_imAdjointVelocityd1);
			this->register_import(m_imAdjointVelocityd2);
			this->register_import(m_imAdjointVelocityd3);
			
			this->register_import(m_imAdjointPressure);
			this->register_import(m_imPressure);
			
			this->register_import(m_imAdjointVelocityGradientd1);
			this->register_import(m_imAdjointVelocityGradientd2);
			this->register_import(m_imAdjointVelocityGradientd3);	
			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		JPrime<TDomain>::
		JPrime(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			

			m_step_length = 1.0;

		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'JPrime'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports
			this->register_import(m_imKinViscosity);			

			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->register_import(m_imAdjointVelocityd1);
			this->register_import(m_imAdjointVelocityd2);
			this->register_import(m_imAdjointVelocityd3);
			
			this->register_import(m_imAdjointPressure);
			this->register_import(m_imPressure);
			
			this->register_import(m_imAdjointVelocityGradientd1);
			this->register_import(m_imAdjointVelocityGradientd2);
			this->register_import(m_imAdjointVelocityGradientd3);
			
			
			this->clear_add_fct();
		}

		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void JPrime<Domain1d>::register_all_funcs(int order,
				int quadOrder)
		{
			//	RegularEdge
			UG_THROW("Not implemented.");
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void JPrime<Domain2d>::register_all_funcs(int order,
				int quadOrder)
		{
			
			register_fe_func<Triangle, DimFEGeometry<dim, 2> >();
			register_fe_func<Quadrilateral, DimFEGeometry<dim, 2> >();
		/* 	if (quadOrder != 2 * order + 1) {
				register_fe_func<Triangle, DimFEGeometry<dim> > ();
				register_fe_func<Quadrilateral, DimFEGeometry<dim> > ();
			} */

		/* 	//	special compiled cases

			//	Triangle
			switch (order)
			{
				case 1:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1> ,
							GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
					register_fe_func<Triangle, FEGeom > (); break;
				}
				case 2:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2> ,
							GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
					register_fe_func<Triangle, FEGeom > (); break;
				}
				case 3:{
					typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3> ,
							GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
					register_fe_func<Triangle, FEGeom > (); break;
				}
				default: register_fe_func<Triangle, DimFEGeometry<dim> > (); break;
			}

			//	Quadrilateral
			switch (order)
			{
				case 1:{
					typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
							ReferenceQuadrilateral, 1> , GaussQuadrature<
							ReferenceQuadrilateral, 3> > FEGeom;
					register_fe_func<Quadrilateral, FEGeom > (); break;
				}
				case 2:{
					typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
							ReferenceQuadrilateral, 2> , GaussQuadrature<
							ReferenceQuadrilateral, 7> > FEGeom;
					register_fe_func<Quadrilateral, FEGeom > (); break;
				}
				case 3:{
					typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<
							ReferenceQuadrilateral, 3> , GaussQuadrature<
							ReferenceQuadrilateral, 11> > FEGeom;
					register_fe_func<Quadrilateral, FEGeom > (); break;
				}
				default: register_fe_func<Quadrilateral, DimFEGeometry<dim> > (); break;
			} */
		}
		#endif
		
		#ifdef UG_DIM_3
		template<>
		void JPrime<Domain3d>::register_all_funcs(int order,
				int quadOrder)
		{
			if (quadOrder != 2 * order + 1) {
				register_fe_func<Tetrahedron, DimFEGeometry<dim> > ();
				register_fe_func<Prism, DimFEGeometry<dim> > ();
				register_fe_func<Pyramid, DimFEGeometry<dim> > ();
				register_fe_func<Hexahedron, DimFEGeometry<dim> > ();
			}

			//	special compiled cases

			//	Tetrahedron
			switch (order) 
			{
				case 1:{
					typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
						ReferenceTetrahedron, 1> , GaussQuadrature<
						ReferenceTetrahedron, 5> > FEGeom;
					register_fe_func<Tetrahedron, FEGeom > (); break;
				}
				case 2:{
					typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
							ReferenceTetrahedron, 2> , GaussQuadrature<
							ReferenceTetrahedron, 5> > FEGeom;
					register_fe_func<Tetrahedron, FEGeom > (); break;
				}
				case 3:{
					typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<
							ReferenceTetrahedron, 3> , GaussQuadrature<
							ReferenceTetrahedron, 7> > FEGeom;
					register_fe_func<Tetrahedron, FEGeom > (); break;
				}
				default: register_fe_func<Tetrahedron, DimFEGeometry<dim> > (); break;

			}

			//	Prism
			switch (order)
			{
				case 1:{
					typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1> ,
							GaussQuadrature<ReferencePrism, 2> > FEGeom;
					register_fe_func<Prism, FEGeom > (); break;
				}
				default:
					register_fe_func<Prism, DimFEGeometry<dim> > (); break;
			}

			//	Pyramid
			switch (order)
			{
				default: register_fe_func<Pyramid, DimFEGeometry<dim> > (); break;
			}

			//	Hexahedron
			switch (order)
			{
				case 1:{
					if (quadOrder == 2){
						typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
						ReferenceHexahedron, 1> , GaussQuadrature<
						ReferenceHexahedron, 2> > FEGeom;
						register_fe_func<Hexahedron, FEGeom > (); break;
					}
					else{
						typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
						ReferenceHexahedron, 1> , GaussQuadrature<
						ReferenceHexahedron, 3> > FEGeom;
						register_fe_func<Hexahedron, FEGeom > (); break;
					}
				}
				case 2:{
					typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
							ReferenceHexahedron, 2> , GaussQuadrature<
							ReferenceHexahedron, 7> > FEGeom;
					register_fe_func<Hexahedron, FEGeom > (); break;
				}
				case 3:{
					typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<
							ReferenceHexahedron, 3> , GaussQuadrature<
							ReferenceHexahedron, 11> > FEGeom;
					register_fe_func<Hexahedron, FEGeom > (); break;
				}
				default: register_fe_func<Hexahedron, DimFEGeometry<dim> > (); break;
			}
		}
		#endif
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void JPrime<TDomain>::register_fe_func()
		{
			ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
			typedef this_type T;
			//static const int refDim = reference_element_traits<TElem>::dim;
			
			this->clear_add_fct(id);

			this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
			this->set_prep_elem_fct(id, &T::template prep_elem<TElem, TFEGeom>);
			this->set_fsh_elem_loop_fct(id, &T::template fsh_elem_loop<TElem, TFEGeom>);

			this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
			this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
			this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
			this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
			this->set_add_rhs_elem_fct(id, &T::template add_rhs_elem<TElem, TFEGeom>);
		}

		////////////////////////////////////////////////////////////////////////////////
		//	explicit template instantiations
		////////////////////////////////////////////////////////////////////////////////

		#ifdef UG_DIM_2
		template class JPrime<Domain2d> ;
		#endif
		#ifdef UG_DIM_3
		template class JPrime<Domain3d> ;
		#endif


}
}

