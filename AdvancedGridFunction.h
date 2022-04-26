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

#ifndef __ADVANCED_GRIDF_FUNCTION
#define __ADVANCED_GRIDF_FUNCTION

#include "lib_disc/function_spaces/grid_function.h"
namespace ug{
  namespace PLaplacian{
    template <typename TDomain, typename TAlgebra>
    class AdvancedGridFunction: public GridFunction <TDomain, TAlgebra>
    {

        public:
            //these functions do the fun
            bool change_storage_type_to_consistent(){
            #ifdef UG_PARALLEL
                this->change_storage_type(PST_CONSISTENT);
            #endif
            return true;
            }
            bool has_storage_type_consistent(){
                return this->has_storage_type(PST_CONSISTENT);
            }
            bool change_storage_type_to_additive(){
            #ifdef UG_PARALLEL
                this->change_storage_type(PST_ADDITIVE);
            #endif
            return true;
            }
            bool has_storage_type_additive(){
                return this->has_storage_type(PST_ADDITIVE);
            }
            bool change_storage_type_to_unique(){
            #ifdef UG_PARALLEL
                this->change_storage_type(PST_UNIQUE);
            #endif
            return true;
            }
            bool has_storage_type_unique(){
                return this->has_storage_type(PST_UNIQUE);
            }

            bool enforce_consistency(){
              change_storage_type_to_unique();
              change_storage_type_to_consistent();
            }

            /// Initializing Constructor
		        AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		        SmartPtr<DoFDistribution> spDoFDistr, bool bManage = true);
         

    	      /// Initializing Constructor using surface dof distribution
		        AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, bool bManage = true);
       

	          /// Initializing Constructor using surface dof distribution on a level
		        AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, int level, bool bManage = true);
       
            /// Initializing Constructor using a grid level
            AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const GridLevel& gl, bool bManage = true);
    };
    /// Initializing Constructor
    template <typename TDomain, typename TAlgebra>
    AdvancedGridFunction<TDomain, TAlgebra>::
		AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
		SmartPtr<DoFDistribution> spDoFDistr, bool bManage): GridFunction<TDomain, TAlgebra>(spApproxSpace, spDoFDistr, bManage)
    {}

    /// Initializing Constructor using surface dof distribution
    template <typename TDomain, typename TAlgebra>
    AdvancedGridFunction<TDomain, TAlgebra>::
		AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, bool bManage ): GridFunction<TDomain, TAlgebra>(spApproxSpace, bManage)
    {}

	  /// Initializing Constructor using surface dof distribution on a level
            template <typename TDomain, typename TAlgebra>
            AdvancedGridFunction<TDomain, TAlgebra>::
		        AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, int level, bool bManage ):GridFunction<TDomain, TAlgebra>(spApproxSpace, level,bManage)
            {}
            /// Initializing Constructor using a grid level
            template <typename TDomain, typename TAlgebra>
            AdvancedGridFunction<TDomain, TAlgebra>::
            AdvancedGridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const GridLevel& gl, bool bManage ):GridFunction<TDomain, TAlgebra>(spApproxSpace, gl,bManage)
            {}
  }
}
#endif
