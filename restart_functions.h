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
#ifndef PLUGIN_PLAPLACIAN_RESTART
#define PLUGIN_PLAPLACIAN_RESTART

#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/normal_calculation_impl.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/disc_util/fv_util.h" //AverageVectors
#include "common/math/math_vector_matrix/math_vector_functions.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

#include "lib_disc/io/vtkoutput.h"

#include <time.h>
namespace ug{
namespace PLaplacian{

///Yes, this is a renamed function...
template <typename TDomain, typename TAlgebra>
static size_t componentName2FctID(const GridFunction<TDomain, TAlgebra>& u,  const char* fctNames, std::vector<size_t> &vFct)
{
	/// translate function names to IDs
		std::vector<std::string> vFctNameStrings;
		std::string fctString(fctNames);
		TokenizeString(fctString, vFctNameStrings, ',');

		 vFct.resize(vFctNameStrings.size());
		for(size_t i = 0; i < vFctNameStrings.size(); ++i)
		{
			vFct[i] = u.fct_id_by_name(vFctNameStrings[i].c_str());
//			std::cout << " Found " << vFctNameStrings[i] << "->" << vFct[i] << std::endl;
		}

		return vFctNameStrings.size();
};

/// Save Nodal positions to a GridFunction, using the domain within the GF
/// fctNames size must be equal to dimensions of the domain
template <typename TDomain, typename TAlgebra>
void SaveNodalPositions2GridFunction(SmartPtr<GridFunction<TDomain, TAlgebra> > inputGF,
								const char* fctNames)
{
    // dimension
	static const int dim = TDomain::dim;

	// translate name to ID
	std::vector<size_t> vFct(dim);
	componentName2FctID(*inputGF, fctNames, vFct);


	// get domain
	SmartPtr<TDomain> domain = inputGF->domain();

	// get dof distribution
	SmartPtr<DoFDistribution> dd = inputGF->dd();

	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	iterator
	typename DoFDistribution::traits<Vertex>::iterator iter, iterEnd;

//	get iterators
	iter = dd->begin<Vertex>(SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<Vertex>(SurfaceView::ALL);

    //	loop all vertices
	for(;iter != iterEnd; ++iter)
	{
		//	get vertex
		Vertex* v = *iter;

		//	index vector
		std::vector<DoFIndex> vMultInd;

		//	compute displacement (must loop all components)
		size_t alpha = 0;
		for(size_t i = 0; i < vFct.size(); ++i)
		{
			//	load indices associated with vertex
			if(dd->inner_dof_indices(v, vFct[i], vMultInd) == 1)
				//myVtxDisp[alpha++] = DoFRef(*inputGF, vMultInd[0]);
                alpha++;
                DoFRef(*inputGF, vMultInd[0]) = aaPos[v][alpha-1];
		}
	}
}//end fnc Nodal2GF


/// Modify domain (vertex positions) by displacement vector 
template <typename TDomain, typename TAlgebra>
void SetDomainCoordinatesFromGF(/*SmartPtr<TDomain> dom,*/
								   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
								   const char* fctNames)
{
        PROFILE_FUNC();
	// dimension
	static const int dim = TDomain::dim;

	// translate name to ID
	std::vector<size_t> vFct(dim);
	componentName2FctID(*u, fctNames, vFct);


	// get domain
	SmartPtr<TDomain> domain = u->domain();

	//UG_ASSERT(dom==domain, "TransformDomainByDisplacement: Domains should match!");

	// get dof distribution
	SmartPtr<DoFDistribution> dd = u->dd();

	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = domain->position_accessor();

//	iterator
	typename DoFDistribution::traits<Vertex>::iterator iter, iterEnd;

//	get iterators
	iter = dd->begin<Vertex>(SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<Vertex>(SurfaceView::ALL);

//	loop all vertices
	for(;iter != iterEnd; ++iter)
	{
		//	get vertex
		Vertex* v = *iter;

		//  extract displacement
		MathVector<dim> myVtxDisp(0.0);

		//	index vector
		std::vector<DoFIndex> vMultInd;

		//	compute displacement (must loop all components)
		size_t alpha = 0;
		for(size_t i = 0; i < vFct.size(); ++i)
		{
			//	load indices associated with vertex
			if(dd->inner_dof_indices(v, vFct[i], vMultInd) == 1)
				myVtxDisp[alpha++] = DoFRef(*u, vMultInd[0]);
			//	UG_THROW("TransformDomainByDisplacement: Need exactly one component, got " <<vMultInd.size() <<"!");
		}
		UG_ASSERT(((alpha ==0) || (alpha==dim)),
				"SetDomainCoordinatesFromGF: No shift or shift with correct dimension!");
		// if (alpha == dim) std::cout << myVtxDisp << std::endl;


		// add displacement to vertex coordinates (on all levels)
		const MultiGrid& mg = *domain->grid();
		while(v){
			//VecAppend(aaPos[v], myVtxDisp);
            VecAssign(aaPos[v],myVtxDisp);
			//GridObject* parent = mg.get_parent(v);
			//if(parent)
			v = dynamic_cast<Vertex*>(mg.get_parent(v));
			// v=0x0, iff parent does not exist or not of vertex-type
		}

	}
}


}//end namespace PLaplacian
}//end namespace ug
#endif