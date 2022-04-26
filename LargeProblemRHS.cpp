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

#include "LargeProblemRHS.h"

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
		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		IdentityDeformationGradient(MathMatrix<dim, dim>& IdGradU, const size_t ip)
		{
			//1.- calculate gradients of each deformation component u1,u2
			MathVector<dim> vdeformation_d1=m_imDeformationVectord1[ip];
			MathVector<dim> vdeformation_d2=m_imDeformationVectord2[ip];
			//2.- assign the gradients to each row using assign(vector, row)
			IdGradU.assign(vdeformation_d1,0);
			IdGradU.assign(vdeformation_d2,1);
			if(this->dim == 3){
				MathVector<dim> vdeformation_d3=m_imDeformationVectord3[ip];
				IdGradU.assign(vdeformation_d3,2);
			}
			//3.- add identity on the diagonal
			for(int i=0; i < dim; ++i)
			{
				IdGradU[i][i] += 1.0;
			}
		}
		//Calculate a matrix using the imports for velocity gradient
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
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
		void LargeProblemRHS<TDomain>::
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
		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		DeformationGradient(MathMatrix<dim, dim>& GradW, const size_t ip)
		{
			//1.- calculate gradients of each deformation component u1,u2
			MathVector<dim> vdeformation_d1=m_imDeformationVectord1[ip];
			MathVector<dim> vdeformation_d2=m_imDeformationVectord2[ip];
			//2.- assign the gradients to each row using assign(vector, row)
			GradW.assign(vdeformation_d1,0);
			GradW.assign(vdeformation_d2,1);
			if(this->dim == 3){
				MathVector<dim> vdeformation_d3=m_imDeformationVectord3[ip];
				GradW.assign(vdeformation_d3,2);
			}
			
		}
		//Calculates the velocities at the ip and puts them into a vector 
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		DeformationVector(MathVector<dim>& w, const size_t ip)
		{
			number w1= m_imDeformationd1[ip];
			number w2= m_imDeformationd2[ip];
			
			w[0]=w1;
			w[1]=w2;
			if(this->dim == 3){
				number w3= m_imDeformationd3[ip];
				w[2]=w3;
			}
		}
		//Calculate a matrix using the imports for adjoint velocity gradient
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
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
		void LargeProblemRHS<TDomain>::
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
		void LargeProblemRHS<TDomain>::
		set_pressure(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imPressure.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_pressure(number val)
		{
			set_pressure(make_sp(new ConstUserNumber<dim>(val)));
		}//***************IMPORT: ADJOINT PRESSURE SCALAR VALUE
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_pressure(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointPressure.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_pressure(number val)
		{
			set_adjoint_pressure(make_sp(new ConstUserNumber<dim>(val)));
		}

		//***************IMPORT: ADJOINT VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d1(number val)
		{
			set_adjoint_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d2(number val)
		{
			set_adjoint_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAdjointVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_d3(number val)
		{
			set_adjoint_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: ADJOINT VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d1(number val)
		{
			set_adjoint_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d2(number val)
		{
			set_adjoint_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imAdjointVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_adjoint_velocity_vector_d3(number val)
		{
			set_adjoint_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		//***************IMPORT: VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d1(number val)
		{
			set_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d2(number val)
		{
			set_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_d3(number val)
		{
			set_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d1(number val)
		{
			set_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d2(number val)
		{
			set_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_velocity_vector_d3(number val)
		{
			set_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		//***************IMPORT: DEFORMATION VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d1(number val)
		{
			set_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d2(number val)
		{
			set_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_d3(number val)
		{
			set_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: DEFORMATION GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		
		//For the d2 
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3 
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}
		
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_kinematic_viscosity(LuaFunctionHandle fct)
		{
			set_kinematic_viscosity(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif

		

		////////////////////////////////////////////////////////////////////////////////
		//	General
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void LargeProblemRHS<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("LargeProblemRHS: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("LargeProblemRHS: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

			//	update assemble functions
			set_assemble_funcs();
		}
		template <typename TDomain>
		void
		LargeProblemRHS<TDomain>::
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
		void LargeProblemRHS<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{

			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			//	prepare geometry for type and order
			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}UG_CATCH_THROW("LargeProblemRHS::prep_elem_loop:"
							" Cannot update Finite Element Geometry.");

			m_imKinViscosity.set_local_ips(geo.local_ips(), geo.num_ip(),false);

			
			m_imDeformationd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imDeformationVectord1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationVectord2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationVectord3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
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
		void LargeProblemRHS<TDomain>::fsh_elem_loop()
		{}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void LargeProblemRHS<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("LargeProblemRHS::prep_elem:"
							" Cannot update Finite Element Geometry.");
			m_imKinViscosity.set_global_ips(geo.global_ips(), geo.num_ip());

			
			m_imDeformationVectord1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationVectord2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationVectord3.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imDeformationd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationd3.set_global_ips(geo.global_ips(), geo.num_ip());
			
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
		void LargeProblemRHS<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	
			//	request geometry
		}
		//  assemble right-hand-side d(i,sh)
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void LargeProblemRHS<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		 
			const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			
			
			number laplacian_power = (m_p-2.0)/2.0;
			
			const number laplacian_power_1 = m_p-1.0;
			const number control_variable = m_control;// pow(1.0/m_control, laplacian_power_1);

			for (size_t ip = 0; ip < geo.num_ip(); ++ip){
								
				
				//Vectors
				MathVector<dim> vVelocity;VecSet(vVelocity,0.0);
				MathVector<dim> vAdjVelocity;VecSet(vAdjVelocity,0.0);
				
				AdjointVelocityVector(vAdjVelocity,ip);
				VelocityVector(vVelocity,ip);
				
				//Declaration of import data structures
				MathMatrix<dim, dim> DefGradient; MatSet(DefGradient, 0.0);
				MathMatrix<dim, dim> IdDeformationGradient;MatSet(IdDeformationGradient,0.0);
				MathMatrix<dim, dim> InverseIdDeformationGradient;MatSet(InverseIdDeformationGradient,0.0);
				number dDF=0.0;
				MathMatrix<dim, dim> AdjVelGrad;MatSet(AdjVelGrad,0.0);
				MathMatrix<dim, dim> VelGrad;MatSet(VelGrad,0.0);
				const number pressure=m_imPressure[ip];
				const number adjoint_pressure=m_imAdjointPressure[ip];
				const number viscosity = m_imKinViscosity[ip];				
				//Calculation
				IdentityDeformationGradient(IdDeformationGradient,ip);
				VelocityGradient(VelGrad,ip);
				AdjointVelocityGradient(AdjVelGrad,ip);
				Inverse(InverseIdDeformationGradient,IdDeformationGradient);
				dDF=Determinant(IdDeformationGradient);
				DeformationGradient(DefGradient, ip);
				
				MathVector<dim> vDeformation;VecSet(vDeformation,0.0);
				DeformationVector(vDeformation,ip);
				
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
						/*				
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
						
						*/
						//Geometrical Terms
						//***************************************************************************
						//******Derivative dV(u,lambda_vol)(delta_u), as part of RHS*****************
						//***************************************************************************
						MathMatrix<dim, dim> TransformedDdelta_u;MatSet(TransformedDdelta_u,0.0);
						MatMultiply(TransformedDdelta_u,InverseIdDeformationGradient,Ddelta_u);
						TrDefInvDdelta_u=Trace(TransformedDdelta_u);
						
						d(d1,sh1) -= m_lambda_vol * dDF *TrDefInvDdelta_u *geo.weight(ip);
						
						//TODO:rmeovve this
						//d(d1,sh1) += TrDdelta_u * geo.weight(ip);
						
						//***************************************************************************
						//******Derivative db(u,lambda_bary)(delta_u), as part of RHS****************
						//***************************************************************************
												
						//Global position of integration point
						const MathVector<dim> vIPPosition=geo.global_ip(ip);
						//Current position
						MathVector<dim> vCurrentPosition;VecSet(vCurrentPosition,0.0);
						VecAdd(vCurrentPosition,vIPPosition,vDeformation);//This is (x+u)

						
						//Barycenter constraints
						//Create shape function vector
						MathVector<dim> vShapeFuncs;VecSet(vShapeFuncs,0.0);
						vShapeFuncs[d1] = geo.shape(ip,sh1);
						MathVector<dim> vShapeFuncsScaled;VecSet(vShapeFuncsScaled,0.0);
						VecScale(vShapeFuncsScaled,vShapeFuncs,dDF);
						MathVector<dim> vCurrentPositionScaled;VecSet(vCurrentPositionScaled,0.0);
						VecScale(vCurrentPositionScaled,vCurrentPosition,TrDefInvDdelta_u*dDF);
						//Create barycenter constraint 1-term
						//addition
						MathVector<dim> vAddition;VecSet(vAddition,0.0);
						VecAdd(vAddition,vShapeFuncsScaled,vCurrentPositionScaled);
						//d(d1,sh1)  += VecDot(m_vLambda_bar , vAddition) * geo.weight(ip);//This performs lambda_var*(integration terms)
						//*********as 2 separate terms
						d(d1,sh1)  -= VecDot(m_vLambda_bar, vShapeFuncs) * dDF * geo.weight(ip);
						//d(d1,sh1) -= m_vLambda_bar[d1]*geo.shape(ip,sh1)* dDF * geo.weight(ip);
						d(d1,sh1)  -= VecDot(m_vLambda_bar, vCurrentPosition) * dDF * TrDefInvDdelta_u * geo.weight(ip);
						//d(d1,sh1)  -= m_vLambda_bar[d1]*vCurrentPosition[d1] * dDF * TrDefInvDdelta_u * geo.weight(ip);
						//***************************************************************************
						//******Derivative dJ(u)(delta_u) p-laplacian, as part of RHS****************
						//***************************************************************************
						number deformation_contraction = MatContraction(DefGradient, Ddelta_u);
						number K1 = MatContraction(DefGradient,DefGradient);// +  m_epsilon*heaviside(4.0 - m_p);
						 
						 
						//pow(base, exponent)
						//Derivative of J wrt u direction delta_u discretized

						//TODO:check this term
						//d(d1,sh1) += pow(K1,laplacian_power)* deformation_contraction* geo.weight(ip);
						for(int d2=0;d2<dim;d2++){
							d(d1, sh1) -= control_variable*std::pow(K1, laplacian_power) * DefGradient(d1, d2) * geo.global_grad(ip, sh1)[d2] * geo.weight(ip) ;
						}
					}
				}
				
				//***************************************************************************
				//******Derivative db(u,lambda_b)(delta_u,delta_lambda), as part of RHS******
				//***************************************************************************
				MathVector<dim> vBarycenterMultipliers;VecSet(vBarycenterMultipliers,0.0);
				vBarycenterMultipliers[0]=m_multiplier_bx;
				vBarycenterMultipliers[1]=m_multiplier_by;
				if (this-> dim == 3){
					vBarycenterMultipliers[2]=m_multiplier_bz;
				}
				
				for (size_t BIndex=0; BIndex < dim ; ++ BIndex){//for recreating the x and y indeces for the vectors Bx and By
					for (int d1 =0; d1 < dim ; ++d1){	
						for (size_t sh1 = 0; sh1 < geo.num_sh(); ++sh1){
							
							number TrTransformedDdeltau=0.0;
							
							const MathVector<dim> vIPPosition=geo.global_ip(ip);
							//Current position, summation of vDeformation + vIPPosition
							MathVector<dim> vCurrentPosition;VecSet(vCurrentPosition,0.0);
							VecAdd(vCurrentPosition,vIPPosition,vDeformation);
							
							MathMatrix<dim, dim> ShapeFuncGradient;MatSet(ShapeFuncGradient,0.0);
							ShapeFuncGradient.assign(geo.global_grad(ip,sh1),d1);
							//ShapeFuncGradient.assign(geo.global_grad(ip,sh1),1);						
							MathMatrix<dim, dim> TransformedDdelta_u;MatSet(TransformedDdelta_u,0.0);
							MatMultiply(TransformedDdelta_u,InverseIdDeformationGradient,ShapeFuncGradient);
							TrTransformedDdeltau=Trace(TransformedDdelta_u);
									
							//Derivative of b(u,lambda_bar)(delta_u,delta_lambda_bar) wrt u and lambda_bar			
							//d(m_index,sh1) += m_multiplier * geo.shape(ip,sh1)*geo.weight(ip)*dDF;//<test_lb,delta_u>*detDF
							d(d1,sh1) -=  vBarycenterMultipliers[BIndex] * vCurrentPosition[BIndex]*TrTransformedDdeltau* geo.weight(ip)*dDF; //<Test_lb,(x+u)>Tr(InvDF*Ddu)*detDF

						}//end sh for
					}
				
					for (size_t sh1 = 0; sh1 < geo.num_sh(); ++sh1){
		
						d(BIndex, sh1) -= vBarycenterMultipliers[BIndex] * geo.shape(ip,sh1) * geo.weight(ip) * dDF;//<test_lb,delta_u>*detDF
					
					}//end sh for				
				}//end index loop, the assembly is done only on the elements of d(index, sh1) that correspond
				//***************************************************************************
				//******Derivative dv(u,lambda_v)(delta_u,delta_lambda), as part of RHS******
				//***************************************************************************
				for (int d1 = 0; d1 < dim; ++d1){
					for (size_t sh1 = 0; sh1 < geo.num_sh(); ++sh1){
					
						number TrDefInvDdelta_u=0.0;
						
						MathMatrix<dim, dim> Ddelta_u;MatSet(Ddelta_u,0.0);
						Ddelta_u.assign(geo.global_grad(ip,sh1),d1);
						
						MathMatrix<dim, dim> TransformedDdelta_u;MatSet(TransformedDdelta_u,0.0);
						MatMultiply(TransformedDdelta_u,InverseIdDeformationGradient,Ddelta_u);
						TrDefInvDdelta_u=Trace(TransformedDdelta_u);

						//Derivative of V(u,lambda_vol)(delta_u,delta_lambda_vol) wrt u and lambda_vol
						//d = Tr(DF-1*Ddelta_u) *detDF*dx
						d(d1,sh1) -= m_multiplier_vol * TrDefInvDdelta_u * dDF*geo.weight(ip);
												
					}//end sh for
				}//end d1 for
				
			}// end for ip
			
	
			
		}
		//  assemble stiffness defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void LargeProblemRHS<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}//end add defect
		
		//	assemble mass jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void LargeProblemRHS<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//  assemble mass-defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void LargeProblemRHS<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}

		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		LargeProblemRHS<TDomain>::
		LargeProblemRHS(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			m_stabilization_scale=1.0;
			
			m_multiplier_vol = 1.0;
			m_multiplier_bx = 1.0;
			m_multiplier_by = 1.0;
			m_multiplier_bz = 1.0;
			
			m_step_length = 1.0;
			m_control = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'LargeProblemRHS'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports
			this->register_import(m_imKinViscosity);			
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imDeformationd1);
			this->register_import(m_imDeformationd2);
			this->register_import(m_imDeformationd3);
			
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
		LargeProblemRHS<TDomain>::
		LargeProblemRHS(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			m_multiplier_vol = 1.0;
			m_multiplier_bx = 1.0;
			m_multiplier_by = 1.0;
			m_multiplier_bz = 1.0;
			m_step_length = 1.0;
			m_control = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'LargeProblemRHS'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports
			this->register_import(m_imKinViscosity);			
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imDeformationd1);
			this->register_import(m_imDeformationd2);
			this->register_import(m_imDeformationd3);
			
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
		void LargeProblemRHS<Domain1d>::register_all_funcs(int order,
				int quadOrder)
		{
			//	RegularEdge
			UG_THROW("Not implemented.");
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void LargeProblemRHS<Domain2d>::register_all_funcs(int order,
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
		void LargeProblemRHS<Domain3d>::register_all_funcs(int order,
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
						ReferenceTetrahedron, 1> > FEGeom;
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
		void LargeProblemRHS<TDomain>::register_fe_func()
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
		template class LargeProblemRHS<Domain2d> ;
		#endif
		#ifdef UG_DIM_3
		template class LargeProblemRHS<Domain3d> ;
		#endif


}
}

