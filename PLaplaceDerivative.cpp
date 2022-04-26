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

#include "PLaplaceDerivative.h"

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
		void PLaplaceDerivative<TDomain>::
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
		//Calculation of the Identity Deformation Gradient
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
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
		void PLaplaceDerivative<TDomain>::
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
		//***************IMPORT: DEFORMATION VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d1(number val)
		{
			set_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d2(number val)
		{
			set_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_d3(number val)
		{
			set_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: DEFORMATION GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		
		//For the d2 
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3 
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}

		////////////////////////////////////////////////////////////////////////////////
		//	General
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void PLaplaceDerivative<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("PLaplaceDerivative: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("PLaplaceDerivative: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

			//	update assemble functions
			set_assemble_funcs();
		}
		template <typename TDomain>
		void
		PLaplaceDerivative<TDomain>::
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
		void PLaplaceDerivative<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{

			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			//	prepare geometry for type and order
			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}UG_CATCH_THROW("PLaplaceDerivative::prep_elem_loop:"
							" Cannot update Finite Element Geometry.");

			m_imDeformationd1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationd2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationd3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			
			m_imDeformationVectord1.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationVectord2.set_local_ips(geo.local_ips(), geo.num_ip(),false);
			m_imDeformationVectord3.set_local_ips(geo.local_ips(), geo.num_ip(),false);
					
		}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::fsh_elem_loop()
		{}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("PLaplaceDerivative::prep_elem:"
							" Cannot update Finite Element Geometry.");
			m_imDeformationVectord1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationVectord2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationVectord3.set_global_ips(geo.global_ips(), geo.num_ip());
			
			m_imDeformationd1.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationd2.set_global_ips(geo.global_ips(), geo.num_ip());
			m_imDeformationd3.set_global_ips(geo.global_ips(), geo.num_ip());

		}
		//  assemble stiffness jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	
			//	request geometry
			const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
			
			for (size_t ip = 0; ip < geo.num_ip(); ++ip){

				MathVector<dim> vDeformation; VecSet(vDeformation, 0.0);	
				MathMatrix<dim, dim> DefGradient; MatSet(DefGradient, 0.0);

				MathMatrix<dim, dim> IdDeformationGradient;MatSet(IdDeformationGradient,0.0);
				MathMatrix<dim, dim> InverseIdDeformationGradient;MatSet(InverseIdDeformationGradient,0.0);
				MathMatrix<dim, dim> NegativeInverseIdDeformationGradient;MatSet(NegativeInverseIdDeformationGradient,0.0);
				number dDF=0.0;
				
				//Calculation
				DeformationVector(vDeformation,ip);
				DeformationGradient(DefGradient, ip);
				number DuContraction=0.0; DuContraction=MatContraction(DefGradient, DefGradient);
				//Add small number to the contraction to avoid singularity
				
				number DuContractionEllipticPart = DuContraction+m_elliptic_epsilon;//*heaviside(m_epsilon-DuContraction);//*heaviside(4.0 - m_p);
				number DuContractionNonLinearPart = DuContraction + m_epsilon*heaviside(4.0 - m_p);//*heaviside(m_epsilon-DuContraction)
				//substraction_tolerance(number p, number cst) for guaranteeing binary zero
				//const number laplacian_power_4=(m_p-4.0)/2.0;
				//const number laplacian_power_2=(m_p-2.0)/2.0;
				//const number laplacian_power_1=(m_p-1.0);
				//With binary zero control
				const number laplacian_power_4=(m_p-4.0)/2.0;
				const number laplacian_power_2=(m_p-2.0)/2.0;
				const number laplacian_power_1=(m_p-1.0);

				//const number K_nonlinear= (m_p-2.0)*pow(DuContractionNonLinearPart, laplacian_power_4);//nonlinear part
				number K_nonlinear= (m_p-2.0)*std::pow(DuContractionNonLinearPart, laplacian_power_4);//nonlinear part
				if (std::abs(m_p-2.0) < 1e-5) K_nonlinear = 0.0;
				number K_elliptic= std::pow(DuContractionEllipticPart, laplacian_power_2);//elliptic par
				if (std::abs(m_p-2.0) < 1e-5) K_elliptic = 1.0;
				const number control_variable = m_control;//pow(1.0/m_control, laplacian_power_1);
							
				IdentityDeformationGradient(IdDeformationGradient,ip);
				Inverse(InverseIdDeformationGradient,IdDeformationGradient);
				MatMultiply(NegativeInverseIdDeformationGradient, InverseIdDeformationGradient, -1.0);
				dDF=Determinant(IdDeformationGradient);
				 

				//TODO:all positives for now
				for(int d1=0; d1 <  dim; ++d1){
					for(size_t sh1=0;sh1 < geo.num_sh();++sh1){
						for(size_t sh2=0; sh2 < geo.num_sh();++sh2){
							for(int d2=0 ; d2< dim ;++d2){
								
								//Second derivative wrt u of J
								MathMatrix<dim,dim> NuGradient;MatSet(NuGradient,0.0);//This coressponds to Dmu_u
								MathMatrix<dim,dim> DeltaGradient;MatSet(DeltaGradient,0.0);//This corresponds to Ddelta_u
								
								NuGradient.assign(geo.global_grad(ip,sh2),d2);								
								DeltaGradient.assign(geo.global_grad(ip,sh1),d1);
								
								//J(d1,sh1,d2,sh2) += control_variable* K_nonlinear* MatContraction(NuGradient, DefGradient)* MatContraction(DeltaGradient, DefGradient) * geo.weight(ip);
								//J(d1,sh1,d2,sh2) += control_variable* K_elliptic* MatContraction(NuGradient, DeltaGradient)*geo.weight(ip);
								J(d1,sh1,d2,sh2) += K_nonlinear* MatContraction(NuGradient, DefGradient)* MatContraction(DeltaGradient, DefGradient) * geo.weight(ip);
								J(d1,sh1,d2,sh2) += K_elliptic*MatContraction(NuGradient, DeltaGradient)*geo.weight(ip);
								
								
								//Second derivative wrt u of V(u,lambda_v)
								number TrDnuDdelta=0.0;
								number TrDnu=0.0;
								number TrDdelta=0.0;
								MathMatrix<dim, dim> TransformedNuGradient;MatSet(TransformedNuGradient,0.0);
								MatMultiply(TransformedNuGradient,InverseIdDeformationGradient,NuGradient);
								
								MathMatrix<dim, dim> NegativeTransformedNuGradient;MatSet(NegativeTransformedNuGradient,0.0);
								MatMultiply(NegativeTransformedNuGradient,NegativeInverseIdDeformationGradient,NuGradient);
								
								MathMatrix<dim, dim> TransformedDeltaGradient; MatSet(TransformedDeltaGradient, 0.0);
								MatMultiply(TransformedDeltaGradient, InverseIdDeformationGradient, DeltaGradient);
								
								MathMatrix<dim, dim> DnuDdelta;MatSet(DnuDdelta, 0.0);
								MatMultiply(DnuDdelta, NegativeTransformedNuGradient, TransformedDeltaGradient);
								
								TrDnuDdelta=Trace(DnuDdelta);//Big trace Tr(InvDf*Dnu*InvDF*Ddelta)
								TrDnu=Trace(TransformedNuGradient);//Tr(invDF*Dnu)
								TrDdelta=Trace(TransformedDeltaGradient);//Tr(invDF*Ddelta)
								
								J(d1,sh1,d2,sh2) += control_variable*m_lambda_vol * TrDnuDdelta * dDF * geo.weight(ip);
								J(d1,sh1,d2,sh2) += control_variable*m_lambda_vol * TrDnu*TrDdelta * dDF * geo.weight(ip);
								
								
								//Second derivative wrt u of b(u,lambda_b)
								
								//Get integration point position, and use interpolated deformation vector u to perform: x + u
								//Global position of integration point
								const MathVector<dim> vIPPosition=geo.global_ip(ip);
								//Current position, summation of vDeformation + vIPPosition
								MathVector<dim> vCurrentPosition;VecSet(vCurrentPosition,0.0);
								VecAdd(vCurrentPosition,vIPPosition,vDeformation);
								
								//Create vectors of shape functions
								//recall Dnu(sh2,d2), as NuGradient
								//recall Ddelta(sh1,d1), as DeltaGradient
								MathVector<dim> vDelta_u;VecSet(vDelta_u,0.0);
								vDelta_u[d1]=geo.shape(ip,sh1);
								
								MathVector<dim> vNu_u;VecSet(vNu_u,0.0);
								vNu_u[d2]=geo.shape(ip,sh2);
								
								J(d1,sh1,d2,sh2) += control_variable*VecDot(m_vLambda_bar,vDelta_u)* TrDnu * dDF * geo.weight(ip);//<vbar,delta_u>*Tr(Inv(DF)Dnu)dDF
								J(d1,sh1,d2,sh2) += control_variable*VecDot(m_vLambda_bar,vNu_u)   * TrDdelta * dDF * geo.weight(ip);//<vbar,mu>*Tr(Inv(DF)*Ddeltau)dDF
								
								number LambdaBarDotXplusU=VecDot(m_vLambda_bar,vCurrentPosition);
								
								J(d1,sh1,d2,sh2) += control_variable*LambdaBarDotXplusU * TrDnu*TrDdelta * dDF * geo.weight(ip);//<lambda_bar,x+u>*Tr(Inv(DF)Dnu)*Tr(Inv(DF)*Ddeltau)dDF
								J(d1,sh1,d2,sh2) += control_variable*LambdaBarDotXplusU * TrDnuDdelta* dDF * geo.weight(ip);//<lambda_bar,x+u>*Tr(-Inv(DF)Dnu*Inv(DF)*Ddeltau)dDF

							}//d2
						}//end sh2
					}//end sh1
				}//d1

			}//end ip
		}
		//  assemble right-hand-side d(i,sh)
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//  assemble stiffness defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}//end add defect
		
		//	assemble mass jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//  assemble mass-defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void PLaplaceDerivative<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}

		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		PLaplaceDerivative<TDomain>::
		PLaplaceDerivative(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			m_stabilization_scale=1.0;
			m_control = 1.0;
			m_step_length = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'PLaplaceDerivative'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports

			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imDeformationd1);
			this->register_import(m_imDeformationd2);
			this->register_import(m_imDeformationd3);
			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		PLaplaceDerivative<TDomain>::
		PLaplaceDerivative(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			m_control = 1.0;
			m_step_length = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'PLaplaceDerivative'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		//	register imports
		
			//For the deformation vectors d1 d2 d3
			this->register_import(m_imDeformationVectord1);
			this->register_import(m_imDeformationVectord2);
			this->register_import(m_imDeformationVectord3);
			
			this->register_import(m_imDeformationd1);
			this->register_import(m_imDeformationd2);
			this->register_import(m_imDeformationd3);
			
			this->clear_add_fct();
		}

		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void PLaplaceDerivative<Domain1d>::register_all_funcs(int order,
				int quadOrder)
		{
			//	RegularEdge
			UG_THROW("Not implemented.");
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void PLaplaceDerivative<Domain2d>::register_all_funcs(int order,
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
		void PLaplaceDerivative<Domain3d>::register_all_funcs(int order,
				int quadOrder)
		{
			typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS< ReferenceTetrahedron, 1> , GaussQuadrature<ReferenceTetrahedron, 8> > HighOrderFEGeom;
			register_fe_func<Tetrahedron, HighOrderFEGeom > ();

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
		void PLaplaceDerivative<TDomain>::register_fe_func()
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
		template class PLaplaceDerivative<Domain2d> ;
		#endif
		#ifdef UG_DIM_3
		template class PLaplaceDerivative<Domain3d> ;
		#endif


}
}

