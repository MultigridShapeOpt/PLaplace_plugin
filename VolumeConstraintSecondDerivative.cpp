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

#include "VolumeConstraintSecondDerivative.h"

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
		void VolumeConstraintSecondDerivative<TDomain>::
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
		
		//Calculation of the  Deformation Gradient
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
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
		void VolumeConstraintSecondDerivative<TDomain>::
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

		////////////////////////////////////////////////////////////////////////////////
		//	IMPORT SETTING FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////

		//***************IMPORT: DEFORMATION VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd1.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d1(number val)
		{
			set_deformation_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd2.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d2(number val)
		{
			set_deformation_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imDeformationd3.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_d3(number val)
		{
			set_deformation_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: DEFORMATION GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord1.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d1(number val)
		{
			set_deformation_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		
		//For the d2 
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord2.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d2(number val)
		{
			set_deformation_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3 
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imDeformationVectord3.set_data(user);
		}
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_deformation_vector_d3(number val)
		{
			set_deformation_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}


		////////////////////////////////////////////////////////////////////////////////
		//	General
		////////////////////////////////////////////////////////////////////////////////
		
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void VolumeConstraintSecondDerivative<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
		//	check number
			if(vLfeID.size() != (size_t)dim)
				UG_THROW("VolumeConstraintSecondDerivative: Needs exactly "<<dim<<" functions.");

			//	check & set order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i].order() < 1)
					UG_THROW("VolumeConstraintSecondDerivative: Adaptive order not implemented.");

		//	remember lfeID;
			m_lfeID = vLfeID[0];

		//	set order
			m_order = vLfeID[0].order();

			//	update assemble functions
			set_assemble_funcs();
		}
		template <typename TDomain>
		void
		VolumeConstraintSecondDerivative<TDomain>::
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
		void VolumeConstraintSecondDerivative<TDomain>::
		prep_elem_loop(const ReferenceObjectID roid, const int si)
		{

			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			//	prepare geometry for type and order
			try{
				geo.update_local(roid, m_lfeID, m_quadOrder);
			}UG_CATCH_THROW("VolumeConstraintSecondDerivative::prep_elem_loop:"
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
		void VolumeConstraintSecondDerivative<TDomain>::fsh_elem_loop()
		{}

		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void VolumeConstraintSecondDerivative<TDomain>::
		prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

			try{
				geo.update(elem, vCornerCoords, m_lfeID, m_quadOrder);
			}
			UG_CATCH_THROW("VolumeConstraintSecondDerivative::prep_elem:"
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
		void VolumeConstraintSecondDerivative<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	

		}
		//  assemble right-hand-side d(i,sh)
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void VolumeConstraintSecondDerivative<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	
		///*
			const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
									
			for (size_t ip = 0; ip < geo.num_ip(); ++ip){
				
				//Declaration of import data structures
				MathMatrix<dim, dim> IdDeformationGradient;MatSet(IdDeformationGradient,0.0);
				MathMatrix<dim, dim> InverseIdDeformationGradient;MatSet(InverseIdDeformationGradient,0.0);
				number dDF=0.0;
			
				//Calculation
				IdentityDeformationGradient(IdDeformationGradient,ip);
				Inverse(InverseIdDeformationGradient,IdDeformationGradient);
				dDF=Determinant(IdDeformationGradient);
				
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
						d(d1,sh1) += m_multiplier * TrDefInvDdelta_u * dDF*geo.weight(ip);
												
					}//end sh for
				}//end d1 for		
			}// end for ip 
		//*/
		}
		//  assemble stiffness defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void VolumeConstraintSecondDerivative<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}//end add defect
		
		//	assemble mass jacobian
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void VolumeConstraintSecondDerivative<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//  assemble mass-defect
		template<typename TDomain>
		template<typename TElem, typename TFEGeom>
		void VolumeConstraintSecondDerivative<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u,
				GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}

		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		VolumeConstraintSecondDerivative<TDomain>::
		VolumeConstraintSecondDerivative(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			m_multiplier = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'VolumeConstraintSecondDerivative'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		
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
		VolumeConstraintSecondDerivative<TDomain>::
		VolumeConstraintSecondDerivative(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_order = 1;
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
			m_multiplier = 1.0;
		//	check number of functions
			if(this->num_fct() != dim)
				UG_THROW("Wrong number of functions: The ElemDisc 'VolumeConstraintSecondDerivative'"
					   " needs exactly "<<dim<<" symbolic function.");
			//Default value assigned	   
		
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
		void VolumeConstraintSecondDerivative<Domain1d>::register_all_funcs(int order,
				int quadOrder)
		{
			//	RegularEdge
			UG_THROW("Not implemented.");
		}
		#endif
		#ifdef UG_DIM_2
		template<>
		void VolumeConstraintSecondDerivative<Domain2d>::register_all_funcs(int order,
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
		void VolumeConstraintSecondDerivative<Domain3d>::register_all_funcs(int order,
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
		void VolumeConstraintSecondDerivative<TDomain>::register_fe_func()
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
		template class VolumeConstraintSecondDerivative<Domain2d> ;
		#endif
		#ifdef UG_DIM_3
		template class VolumeConstraintSecondDerivative<Domain3d> ;
		#endif


}
}

