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
 
 // for various user data
 //also taken from restricted_deformation_elasticity.cpp
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

#include "NavierStokes.h"
 
 
 namespace ug{
	 namespace PLaplacian{

		////////////////////////////////////////////////////////////////////////////////
		//	general
		////////////////////////////////////////////////////////////////////////////////
		template<typename TDomain>
		void NavierStokes<TDomain>::set_quad_order(size_t order)
		{
			m_quadOrder = order;
			m_bQuadOrderUserDef = true;
		}
		
		template<typename TDomain>
		void NavierStokes<TDomain>::
		prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
		{
			if(bNonRegularGrid)
				UG_THROW("NavierStokes: only implemented for regular grids.");

			//	check number
			if(vLfeID.size() != dim+1)
				UG_THROW("NavierStokes: Needs exactly "<<dim+1<<" functions.");

			for(int d = 1; d < dim; ++d)
				if(vLfeID[0] != vLfeID[d])
					UG_THROW("NavierStokes: trial spaces for velocity expected to be"
					" identical for all velocity components.");

			//	remember lfeID;
			m_vLFEID = vLfeID[0];
			m_pLFEID = vLfeID[dim];

			if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_vLFEID.order()+1;

			//	update assemble functions
			register_all_funcs(m_vLFEID, m_pLFEID, m_quadOrder);
		}
		// FUNCTIONS THAT NEED REGISTERING
		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_nonlinear(bool option)
		{
			this->m_bNonLinear=option;
		}
		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_picard(bool option)
		{
			this->m_bPicard=option;
		}

		////IMPORTS
		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imKinViscosity.set_data(user);
		}
		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_kinematic_viscosity(number val)
		{
			if(val == 0.0) set_kinematic_viscosity(SmartPtr<CplUserData<number, dim> >());
			else set_kinematic_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}
		#ifdef UG_FOR_LUA
		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_kinematic_viscosity(const char* fctName)
		{
			set_kinematic_viscosity(LuaUserDataFactory<number,dim>::create(fctName));
		}

		template<typename TDomain>
		void NavierStokes<TDomain>::
		set_kinematic_viscosity(LuaFunctionHandle fct)
		{
			set_kinematic_viscosity(make_sp(new LuaUserData<number,dim>(fct)));
		}
		#endif
		
		////////////////////////////////////////////////////////////////////////////////
		// Assembling functions
		////////////////////////////////////////////////////////////////////////////////
		//Element looping functions
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::prep_elem_loop(const ReferenceObjectID roid, const int si)
		{
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try
			{
				vgeo.update_local(roid, m_vLFEID, m_quadOrder);
				pgeo.update_local(roid, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("FluidOptim::NavierStokes: "
			"Cannot update Finite Element Geometry.");
			m_imKinViscosity.set_local_ips(vgeo.local_ips(), vgeo.num_ip(),false);
		}
		
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
		{
			m_pElem = elem;

			// 	Update Geometry for this element
			DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
			try{
				vgeo.update(elem, vCornerCoords, m_vLFEID, m_quadOrder);
				pgeo.update(elem, vCornerCoords, m_pLFEID, m_quadOrder);
			}
			UG_CATCH_THROW("NavierStokes: Cannot update "
							"Finite Element Geometry.");
			static const int refDim = TElem::dim;
			m_imKinViscosity.set_global_ips(vgeo.global_ips(), vgeo.num_ip());

		}
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::fsh_elem_loop()
		{
		}
		///	assembles the local stiffness matrix 
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::
		add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		 //request geometry
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);
										
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
												
			////////////////////////////////////////////////////
			// Diffusive Term (Momentum Equation)
			////////////////////////////////////////////////////

				const number scale = m_imKinViscosity[ip]* vgeo.weight(ip);
						
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
							for (int udim = 0; udim < dim; ++udim) {

							J(vdim, vsh, vdim, ush) +=  scale
										  * vgeo.global_grad(ip, ush)[udim]
										  * vgeo.global_grad(ip, vsh)[udim];
							}
						}
					}
				}


				if(m_bNonLinear == true)
				{
				////////////////////////////////////////////////////
				// Linearization Terms
				////////////////////////////////////////////////////
							
				////////////////////////////////////////////////////
				//Vector-Convection Matrix (N)
				////////////////////////////////////////////////////
				//interpolate velocity at ip, assign quasi-laplacian terms in d1,d1

				// Interpolate Velocity at ip
				MathVector<dim> Vel;
				for(int d = 0; d < dim; ++d){
					Vel[d] = 0.0;
					for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
						Vel[d] += u(d, sh) * vgeo.shape(ip, sh);
				}

				// linearization of u \nabla u w.r.t second u (i.e. keeping first as old iterate)
				for (int vdim = 0; vdim < dim; ++vdim){
					for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
						for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
							for (int udim = 0; udim < dim; ++udim) {
								J(vdim, vsh, vdim, ush) += Vel[udim]
											 * vgeo.global_grad(ip, ush)[udim]
											 * vgeo.shape(ip, vsh)
											 * vgeo.weight(ip);
							}
						}
					}
				}
							
							
					if(m_bPicard == false)
					{
					////////////////////////////////////////////////////
					//Newton-Derivative Matrix 
					////////////////////////////////////////////////////
								
					// Interpolate Functional Matrix of velocity at ip
					MathMatrix<dim, dim> gradVel;
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 <dim; ++d2){
							gradVel(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
								gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
						}
					}

					// linearization of u \nabla u w.r.t first u (i.e. keeping second as old iterate)
					for (int vdim = 0; vdim < dim; ++vdim){
						for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
							for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
								for (int udim = 0; udim < dim; ++udim) {
									J(vdim, vsh, udim, ush) += gradVel(vdim, udim)
												* vgeo.shape(ip, ush)
												* vgeo.shape(ip, vsh)
												* vgeo.weight(ip);
								}
							}
						}
					}
							
					}//end if picard														
				}//if nonlinear

				////////////////////////////////////////////////////
				// Pressure Term (Momentum Equation)
				////////////////////////////////////////////////////
				for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
					for (int vdim = 0; vdim < dim; ++vdim){
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
							J(vdim, vsh, _P_, psh) -= pgeo.shape(ip, psh)
										* vgeo.global_grad(ip, vsh)[vdim]
										* vgeo.weight(ip);
						}
					}
				}
							
				////////////////////////////////////////////////////
				// Continuity Equation (conservation of mass)
				////////////////////////////////////////////////////

				for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
					for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
						for (int udim = 0; udim < dim; ++udim) {
							J(_P_, psh, udim, ush) +=
										  vgeo.global_grad(ip, ush)[udim]
										* pgeo.shape(ip, psh)
										* vgeo.weight(ip);
						}
					}
				}

	
						
						
			}//end ip

			if(m_stab_type ==  0.0 && m_stabParam != 0.0)
			{
			const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());

				for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
					for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
						for (size_t psh2 = 0; psh2 < pgeo.num_sh(); ++psh2){
							J(_P_, psh,_P_, psh2) += scale
										* VecDot(pgeo.global_grad(ip, psh), pgeo.global_grad(ip, psh2))
										* pgeo.weight(ip);
						}
					}
				}
							
						

																		
			}//case pressure_gradient
			else if(m_stab_type != 0.0)
			{
				for(size_t ip=0; ip<pgeo.num_ip();++ip)
				{
					const number scale = m_stab_type/m_imKinViscosity[ip];
									
					for(size_t psh1=0;psh1<pgeo.num_sh();++psh1)
					{	
						for(size_t psh2=0;psh2<pgeo.num_sh();psh2++){
							J(_P_,psh1,_P_,psh2)+= scale*(pgeo.shape(ip,psh1)*pgeo.shape(ip,psh2)-1.0/(pgeo.num_sh()*pgeo.num_sh()))*pgeo.weight(ip);
						}
					}
				}						
			}//case pressure_projection

		}//end add_jac
		
		
		///	assembles the stiffness part of the local defect, for nonlinear stuff
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::
		add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		// request geometry
			const DimFEGeometry<dim>& vgeo = GeomProvider<VGeom>::get(m_vLFEID, m_quadOrder);
			const DimFEGeometry<dim>& pgeo = GeomProvider<PGeom>::get(m_pLFEID, m_quadOrder);

			// loop integration points, note: pgeo and vgeo have same ip
			for (size_t ip = 0; ip < vgeo.num_ip(); ++ip){
					
					// 	Interpolate Functional Matrix of velocity at ip
					MathMatrix<dim, dim> gradVel;MatSet(gradVel,0.0);
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 <dim; ++d2){
							gradVel(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
								gradVel(d1, d2) += u(d1, sh) * vgeo.global_grad(ip, sh)[d2];
						}
					}

					number divu = 0.0;
					for (size_t ush = 0; ush < vgeo.num_sh(); ++ush){
						for (int udim = 0; udim < dim; ++udim) {
							divu += u(udim, ush) * vgeo.global_grad(ip, ush)[udim];
						}
					}
					
					////////////////////////////////////////////////////
					// Diffusive Term (Momentum Equation)
					////////////////////////////////////////////////////
					const number scale = m_imKinViscosity[ip] * vgeo.weight(ip);
					for (int vdim = 0; vdim < dim; ++vdim){
						for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
							for (int udim = 0; udim < dim; ++udim) {
								d(vdim, vsh) +=  scale * gradVel(vdim, udim)
										       * vgeo.global_grad(ip, vsh)[udim];
							}
						}
					}
					
					if(m_bNonLinear)
					{
						// 	Interpolate Velocity at ip
							MathVector<dim> Vel;VecSet(Vel,0.0);
							for(int d1 = 0; d1 < dim; ++d1){
								Vel[d1] = 0.0;
								for(size_t sh = 0; sh < vgeo.num_sh(); ++sh)
									Vel[d1] += u(d1, sh) * vgeo.shape(ip, sh);
							}

							MathVector<dim> convFlux;VecSet(convFlux,0.0);
							MatVecMult(convFlux, gradVel, Vel);

							for (int vdim = 0; vdim < dim; ++vdim){
								for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
									d(vdim, vsh) += convFlux[vdim]
										      * vgeo.shape(ip, vsh) * vgeo.weight(ip);
								}
							}
						
						
					}//end if nonlinear
				
					////////////////////////////////////////////////////
					// Pressure Term (Momentum Equation)
					////////////////////////////////////////////////////
					number pressure = 0.0;
					for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
						pressure += u(_P_, psh) * pgeo.shape(ip, psh);
					}
						for (size_t vsh = 0; vsh < vgeo.num_sh(); ++vsh){
							for (int vdim = 0; vdim < dim; ++vdim){
								d(vdim, vsh) -= pressure
											   * vgeo.global_grad(ip, vsh)[vdim]
											   * vgeo.weight(ip);
							}
						}
					

					////////////////////////////////////////////////////
					// Continuity Equation (conservation of mass)
					////////////////////////////////////////////////////

						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
							d(_P_, psh) += divu * pgeo.shape(ip, psh)* vgeo.weight(ip);
						}
					
				}

				// stabilization
				if(m_stab_type == 0.0 && m_stabParam != 0.0)
				{

					const number scale = m_stabParam * ElementDiameterSq<GridObject, TDomain>(*m_pElem, *this->domain());

						for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){

							MathVector<dim> pressGrad; VecSet(pressGrad, 0.0);
							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
								VecScaleAppend(pressGrad, u(_P_, psh), pgeo.global_grad(ip, psh));

							for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
								d(_P_, psh) += scale
											   * VecDot(pgeo.global_grad(ip, psh), pressGrad)
											   * pgeo.weight(ip);
							}
						}

				}//case pressure_gradient
				else if (m_stab_type != 0.0)
				{
					number pressure_average = 0.0;
					for (size_t psh = 0; psh < pgeo.num_sh(); ++psh) {
						pressure_average +=  1.0/pgeo.num_sh()*u(_P_, psh);
					}
					for (size_t ip = 0; ip < pgeo.num_ip(); ++ip){
						
						const number scale = m_stab_type/m_imKinViscosity[ip];

								
						number pressure = 0.0;
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh)
							pressure += u(_P_, psh) * pgeo.shape(ip, psh);
							
						for (size_t psh = 0; psh < pgeo.num_sh(); ++psh){
							//d(_P_, psh) += m_imKinViscosity[ip]*pgeo.shape(ip,psh)*pressure*(1.0 - 1.0/(pgeo.num_sh()*pgeo.num_sh())) * pgeo.weight(ip)*dDF;
							//d(_P_, psh) += pressure*(pgeo.shape(ip,psh)*pgeo.shape(ip,psh) - (1.0/(pgeo.num_sh()*pgeo.num_sh()))) * pgeo.weight(ip)*dDF;
							d(_P_, psh) += scale*(pressure * pgeo.shape(ip,psh) - pressure_average/pgeo.num_sh()) * pgeo.weight(ip);
						}
					}//ip end
				}//case pressure_projection

		}
		///	assembles the local right hand side
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
			
		}
		//Mass Matrix methods will be assembled if there is time dependency
		///	assembles the local mass matrix using a finite volume scheme
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		}
		///	assembles the mass part of the local defect
		template<typename TDomain>
		template<typename TElem, typename VGeom, typename PGeom>
		void NavierStokes<TDomain>::add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
		{
		}
		
		////////////////////////////////////////////////////////////////////////////////
		//	Constructor
		////////////////////////////////////////////////////////////////////////////////

		template<typename TDomain>
		NavierStokes<TDomain>::
		NavierStokes(const char* functions, const char* subsets)
		 : IElemDisc<TDomain>(functions,subsets)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
		//Default value assigned	   
		m_stabParam=0.0;// no stab
		m_stab_type=0.0;
				
		//	register imports
			this->register_import(m_imKinViscosity);			
			this->clear_add_fct();
		}
		
		template<typename TDomain>
		NavierStokes<TDomain>::
		NavierStokes(const std::vector<std::string>& vFct,
                                    const std::vector<std::string>& vSubset)
		 : IElemDisc<TDomain>(vFct,vSubset)
		{
		//	set defaults
			m_bQuadOrderUserDef = false;
			m_quadOrder = -1;
			
		//	check number of functions
			if(this->num_fct() != dim+1)
				UG_THROW("Wrong number of functions: The ElemDisc 'NavierStokes'"
					   " needs exactly "<<dim+1<<" symbolic function.");
			//Default value assigned	   
			m_stabParam=0.0;
			m_stab_type=0.0;

		//	register imports
			this->register_import(m_imKinViscosity);			
			this->clear_add_fct();
		}

		////////////////////////////////////////////////////////////////////////////////
		//	register assemble functions
		////////////////////////////////////////////////////////////////////////////////
		#ifdef UG_DIM_1
		template<>
		void NavierStokes<Domain1d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			UG_THROW("Not implemented.");
		}
		#endif

		#ifdef UG_DIM_2
		template<>
		void NavierStokes<Domain2d>::
		register_all_funcs(const LFEID& vLfeID, const LFEID& pLfeID, const int quadOrder)
		{
			typedef DimFEGeometry<dim> FVGeom;
			register_func<Triangle, FVGeom, FVGeom >();
			register_func<Quadrilateral, FVGeom, FVGeom >();
		}
		#endif

		#ifdef UG_DIM_3
		template<>
		void NavierStokes<Domain3d>::
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
		void NavierStokes<TDomain>::register_func()
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
		template class NavierStokes<Domain2d>;
		#endif
		#ifdef UG_DIM_3
		template class NavierStokes<Domain3d>;
		#endif
	 }//namespace NavierStokes
 }//namespace ug4
