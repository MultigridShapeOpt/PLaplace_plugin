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
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"

#include "NavierStokesAdjoint.h"

namespace ug{
	namespace PLaplacian{
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		update_geo_elem(TBaseElem* elem, DimFEGeometry<dim>& vgeo, DimFEGeometry<dim>& pgeo)
		{
			
			SmartPtr<TDomain> dom = this->domain();

			typedef typename IElemDisc<TDomain>::domain_type::position_accessor_type
					position_accessor_type;
			const position_accessor_type& aaPos = dom->position_accessor();

			//	coord and vertex array
			MathVector<dim> coCoord[domain_traits<dim>::MaxNumVerticesOfElem];

			//	get vertices and extract corner coordinates
			const size_t numVertices = elem->num_vertices();
			for (size_t i = 0; i < numVertices; ++i) {
				coCoord[i] = aaPos[elem->vertex(i)];
			};

			//	prepare geometry for type and order
		   	try{
				vgeo.update(elem, &(coCoord[0]), m_vLFEID, m_quadOrder);
				pgeo.update(elem, &(coCoord[0]), m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("AdjointNavierStokes::update_geo_elem:"
							" Cannot update Finite Element Geometry.");
		}
		////////////////////////////////////////////////////////////////////////////////
		//NON-REGISTERED FUNCTIONS: THEY USE IMPORTS TO CALCULATE, NOT PUBLIC
		////////////////////////////////////////////////////////////////////////////////
		
		//Calculate a matrix using the imports for velocity gradient
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
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
		void NavierStokesAdjoint<TDomain>::
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
		
		////////////////////////////////////////////////////////////////////////////////
		//	IMPORT SETTING FUNCTIONS
		////////////////////////////////////////////////////////////////////////////////
		
		//Let's not use the LUA handlers for these imports, since we pass a CplUserData derived class directly :)
		
		//***************IMPORT: VELOCITY VECTOR COMPONENTS
		//d1
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d1(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd1.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d1(number val)
		{
			set_velocity_d1(make_sp(new ConstUserNumber<dim>(val)));
		}		
		//d2
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d2(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd2.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d2(number val)
		{
			set_velocity_d2(make_sp(new ConstUserNumber<dim>(val)));
		}
		//d3
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d3(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imVelocityd3.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_d3(number val)
		{
			set_velocity_d3(make_sp(new ConstUserNumber<dim>(val)));
		}
		//***************IMPORT: VELOCITY GRADIENT ROWS
		//For the d1
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d1(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd1.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d1(number val)
		{
			set_velocity_vector_d1(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d2
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d2(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd2.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d2(number val)
		{
			set_velocity_vector_d2(make_sp(new ConstUserVector<dim>(val)));
		}
		//For the d3
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d3(SmartPtr<CplUserData<MathVector<dim>,dim>> user)
		{
			m_imVelocityGradientd3.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_velocity_vector_d3(number val)
		{
			set_velocity_vector_d3(make_sp(new ConstUserVector<dim>(val)));
		}

		////////////////////////////////////////////////////////////////////////////////
		//	general
		////////////////////////////////////////////////////////////////////////////////
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		prepare_setting(const std::
		vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
			if(bNonRegularGrid)
				UG_THROW("NavierStokesAdjoint: only implemented for regular grids.");

			//	check number
			if(vLfeID.size() != dim+1)
				UG_THROW("NavierStokesAdjoint: Needs exactly "<<dim+1<<" functions.");

			for(int d = 1; d < dim; ++d)
				if(vLfeID[0] != vLfeID[d])
					UG_THROW("NavierStokesAdjoint: trial spaces for velocity expected to be"
					" identical for all velocity components.");

			//	remember lfeID;
			m_vLFEID = vLfeID[0];
			m_pLFEID = vLfeID[dim];

			if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_vLFEID.order()+1;

			//	update assemble functions
			register_all_funcs(m_vLFEID, m_pLFEID, m_quadOrder);
		}
/*		
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void NavierStokesAdjoint<TDomain>::
		set_kinematic_viscosity(LuaFunctionHandle fct)
		{
			set_kinematic_viscosity(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif
*/		
		////////////////////////////////////////////////////////////////////////////////
		// Assembling functions
		////////////////////////////////////////////////////////////////////////////////
		
		//Element looping functions
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try
			{
				vgeo.update_local(roid, m_vLFEID, m_quadOrder);
				pgeo.update_local(roid, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("PLaplacian::NavierStokesAdjoint: "
			"Cannot update Finite Element Geometry.");
			//m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);

			
			m_imVelocityGradientd1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityGradientd2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityGradientd3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			
			m_imVelocityd1.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityd2.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
			m_imVelocityd3.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, 
											   const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			m_pElem = elem;

			// 	Update Geometry for this element
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try{
				vgeo.update(elem, vCornerCoords, m_vLFEID, m_quadOrder);
				pgeo.update(elem, vCornerCoords, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("NavierStokesAdjoint: Cannot update "
							"Finite Element Geometry.");

			//m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());

			
			m_imVelocityGradientd1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityGradientd2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityGradientd3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
			m_imVelocityd1.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityd2.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			m_imVelocityd3.set_global_ips(vgeo.global_ips(), vgeo.num_ip());
			
		}
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::fsh_elem_loop()
		{
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			//	request geometry
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
													
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
				
			
				////////////////////////////////////////////////////
				// Stiffness/Diffusion Matrix Terms
				////////////////////////////////////////////////////
				const number scale = m_imKinViscosity* vgeo.weight(ip);
				
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
							for (int udim = 0; udim < dim; ++udim) {
								J(vdim, vsh, vdim, ush) +=  scale
											 * vgeo.global_grad(ip,ush)[udim]
											 * vgeo.global_grad(ip,vsh)[udim];

							}
						}
					}
				}
				
				/////////////////////////////////////////////////////////////
				// Newton-Derivative Matrix, uses Gradient of Velocity at ip 
				/////////////////////////////////////////////////////////////
				// 	Interpolate Functional Matrix of velocity at ip
				MathMatrix<dim, dim> GradVel;MatSet(GradVel,0.0);
				VelocityGradient(GradVel,ip);
																	
				
				for(int vdim=0 ; vdim < dim; vdim++){
					for(size_t vsh=0 ; vsh < vgeo.num_sh() ; vsh++){
						for(size_t ush=0 ; ush<vgeo.num_sh() ; ush++){
							for(int udim=0 ; udim < dim ; udim++){
								J(vdim,vsh,udim,ush) += GradVel(udim,vdim)
											 			*vgeo.shape(ip,ush)
														*vgeo.shape(ip,vsh)
														*vgeo.weight(ip); 
							}
						}
					}

				}//end sh1
				////////////////////////////////////////////////////
				// Vector-Convection Matrix, uses Velocity value at ip
				////////////////////////////////////////////////////
				
				// 	Interpolate Velocity at ip
				MathVector<dim> Vel;VecSet(Vel,0.0);
				VelocityVector(Vel,ip);
				for(int vdim=0; vdim < dim ; vdim++){
					for(size_t vsh=0; vsh < vgeo.num_sh() ;++vsh){
						for(size_t ush=0; ush < vgeo.num_sh(); ush++){
							
							
							MathMatrix<dim,dim> GradientQ;MatSet(GradientQ,0.0);
							GradientQ.assign(vgeo.global_grad(ip,vsh),vdim);
	
							//Perform GradQInvDF*V
							MathVector<dim> vTemp1;VecSet(vTemp1,0.0);//MatVec multiplication
							MatVecMult(vTemp1,GradientQ,Vel);
							//Perform dot product inside loop, 
							for(int udim=0; udim < dim; udim++){

								//Performs (Grad_Test_Q*V).Q
								MathVector<dim> vShapes;VecSet(vShapes,0.0);vShapes[udim]=vgeo.shape(ip,ush);
								number vec_prod = VecProd(vShapes,vTemp1);
								
								J(vdim,vsh,udim,ush) += vec_prod*vgeo.weight(ip);
							} 
 							
							

						}
					}
				}//end sh1
				
				
				////////////////////////////////////////////////////
				// Adjoint Pressure Term H
				////////////////////////////////////////////////////
				
				for(int vdim=0 ; vdim < dim ; ++vdim){
					for(size_t vsh=0;vsh<vgeo.num_sh();vsh++){
						for(size_t psh=0 ; psh < pgeo.num_sh() ; ++psh){
										
							
							//Trace(GradientQ)
							MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
							GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);
							//trace calculation
							number TrAdjVel = Trace(GradientAdjoint);
							//assign to matrix
							J(vdim,vsh,_P_,psh) -= TrAdjVel*pgeo.shape(ip,psh)*vgeo.weight(ip);

							
						}//psh
					}//vdim		
				}//vsh
				
				////////////////////////////////////////////////////
				// Continuity Equation
				////////////////////////////////////////////////////
				//calculate the trace before assigning value
				for(size_t psh=0; psh < pgeo.num_sh() ; ++psh){
					for(size_t vsh=0; vsh < vgeo.num_sh() ; ++vsh){
						for(size_t vdim=0; vdim < (size_t) dim; ++vdim){
								
							//Trace(GradientQ*InvDF)
							MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
							GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);

							number TrAdjVel = Trace(GradientAdjoint);
							
							//assign value to matrix
							J(_P_,psh,vdim,vsh) += TrAdjVel*pgeo.shape(ip,psh)*vgeo.weight(ip);

						}
					}
				}
			}//end ip
			
			////////////////////////////////////////////////////
			// Stabilization Term
			////////////////////////////////////////////////////
			if(m_stab_type == 0.0 && m_stabParam != 0.0)
			{
				const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());
				for(size_t ip=0; ip<pgeo.num_ip();++ip){
									
					for(size_t psh1=0;psh1<pgeo.num_sh();++psh1){
						for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){

							J(_P_,psh1,_P_,psh2)+= scale*VecDot(pgeo.global_grad(ip, psh1), pgeo.global_grad(ip, psh2))*pgeo.weight(ip);
						}
					}
				}//end pressure ip		
			}//end stab param !=0 
			else if(m_stab_type != 0.0)
			{
				for(size_t ip=0; ip<pgeo.num_ip();++ip)
				{		
					const number scale = m_stab_type/m_imKinViscosity;
									
					for(size_t psh1=0;psh1<pgeo.num_sh();++psh1)
					{	
						for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){
							J(_P_,psh1,_P_,psh2)+= scale*(pgeo.shape(ip,psh1)*pgeo.shape(ip,psh2)-1.0/(pgeo.num_sh()*pgeo.num_sh()))*pgeo.weight(ip);						
						}					
					}
				}						
			}//case pressure_projection
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::
		add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
													
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
				//1.-Need Gradient of the velocity field (matrix) and the Deformation Gradient
				MathMatrix<dim, dim> VelGrad;MatSet(VelGrad,0.0);			
				VelocityGradient(VelGrad,ip);
				
				const number scale= m_imKinViscosity*vgeo.weight(ip);
					
				//2.- Assign value to d
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						
						//different approach now
						MathMatrix<dim, dim> GradientAdjoint;MatSet(GradientAdjoint,0.0);
						GradientAdjoint.assign(vgeo.global_grad(ip,vsh),vdim);
						number contraction_VelAdjoint = MatContraction(VelGrad,GradientAdjoint);
						//assign to rhs
						d(vdim, vsh) += scale*contraction_VelAdjoint;

					}
				}
				
			}
		}
		///	assembles the local right hand side
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		//Mass Matrix methods will be assembled if there is time dependency
		///	assembles the local mass matrix using a finite volume scheme
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::
		add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		///	assembles the mass part of the local defect
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::
		add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{	}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		NavierStokesAdjoint<TDomain>::
		NavierStokesAdjoint(const char* functions, const char* subsets)
		:IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			m_stab_type=0.0;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokesAdjoint'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;// no stab
		//	register imports
			//this->register_import(m_imKinViscosity);			

			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		NavierStokesAdjoint<TDomain>::
		NavierStokesAdjoint(const std::vector<std::string>& vFct,
                      const std::vector<std::string>& vSubset)
		:IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults		
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokesAdjoint'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;// no stab
			m_stab_type=0.0;
		//	register imports
			//this->register_import(m_imKinViscosity);
			
			this->register_import(m_imVelocityGradientd1);
			this->register_import(m_imVelocityGradientd2);
			this->register_import(m_imVelocityGradientd3);
			
			this->register_import(m_imVelocityd1);
			this->register_import(m_imVelocityd2);
			this->register_import(m_imVelocityd3);
			
			this->clear_add_fct();
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void NavierStokesAdjoint<Domain1d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			UG_THROW("Not implemented.");
		}
		#endif

		#ifdef UG_DIM_2
		template<>
		void NavierStokesAdjoint<Domain2d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Triangle, FVGeom, FVGeom >();
			register_func<Quadrilateral, FVGeom, FVGeom >();
		}
		#endif

		#ifdef UG_DIM_3
		template<>
		void NavierStokesAdjoint<Domain3d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Tetrahedron, FVGeom, FVGeom >();
			register_func<Prism, FVGeom, FVGeom >();
			register_func<Hexahedron, FVGeom, FVGeom >();
		}
		#endif
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokesAdjoint<TDomain>::register_func()
		{
			ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
			typedef this_type T;

			this->clear_add_fct(id);
			this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, VGeom, PGeom>);
			this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, VGeom, PGeom>);
			this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, VGeom, PGeom>);
			this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, VGeom, PGeom>);
			this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, VGeom, PGeom>);
			this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, VGeom, PGeom>);
			this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, VGeom, PGeom>);
			this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, VGeom, PGeom>);
		}
		////////////////////////////////////////////////////////////////////////////////
		//	explicit template instantiations
		////////////////////////////////////////////////////////////////////////////////

		#ifdef UG_DIM_2
		template class NavierStokesAdjoint<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class NavierStokesAdjoint<Domain3d>;
		#endif
	}//end namespace FluidOpt
}//end namespace ug
