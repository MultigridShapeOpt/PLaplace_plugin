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
#ifndef PLUGIN_PLAPLACIAN_FUNCTIONS
#define PLUGIN_PLAPLACIAN_FUNCTIONS

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



template <typename TGridFunction>
number PNormToP(TGridFunction* spGridFct,
				const char* vCmp,
				const char* IntegrationSubsets,                           
				int quadOrder,
				const number norm_order){

	number scalar_product=0.0;//final result

	static const int dim = TGridFunction::dim;
	static const int WorldDim = TGridFunction::dim;
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object grid_base_object;
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
	
	std::vector<DoFIndex> multInd;
	SubsetGroup intSSGrp(spGridFct->domain()->subset_handler());
	intSSGrp.add(TokenizeString(IntegrationSubsets));
	
	//std::cout<<"L2VecProd: 1.- DEFINITION OF TYPES AND INITIALIZATION DONE\n";
	
	//	Create function group from GridFunction components
	const FunctionGroup vFctID = spGridFct->fct_grp_by_name(vCmp);//should be of size 1 (only 1 function)

	std::vector<LFEID> vLFEID;
	//Extract and store the LFE id from the FunctionPattern through the FunctionGroup
	for(size_t fct = 0; fct < vFctID.size(); ++fct){
		vLFEID.push_back(spGridFct->lfeid(vFctID[fct]));
	}
	

	//std::cout<<"L2VecProd: 2.- CREATION OF FINITE ELEMENT IDS AND EXTRACTION OF COMPONENT DONE\n";
	//std::cout<<"L2VecProd: 3.- START OF LOOP FOR SUBSETS\n";
	//Loop subsets assigned to the integration subset group
	for(size_t i = 0; i < intSSGrp.size(); ++i)
	{
		//Store subset index
		const int si = intSSGrp[i];
		const_iterator iterBegin = spGridFct->template begin<grid_base_object>(si);
		const_iterator iterEnd = spGridFct->template end<grid_base_object>(si);
		const_iterator iter = iterBegin;
		typename domain_traits<TGridFunction::dim>::position_accessor_type& aaPos= spGridFct->domain()->
																				   position_accessor();
		
		typedef typename domain_traits<dim>::element_type Element;

		//vector of MathVectors sized dim, will store the corners of each face
		std::vector<MathVector<WorldDim> > vCorner;
		std::vector<int> vSubsetIndex;
		//For storing a Jacobian per ip for a given element
		//std::vector<MathMatrix<dim, WorldDim> > vJT;
		std::vector<MathMatrix<dim, WorldDim> > vElemJT;

		for(; iter != iterEnd; ++iter)
		{
			//Pointer to element obtained from dereferencing the iterator...see note above
			Element* pElem = *iter;
		
			CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
			
			//This tells us which kind of geometric object it is, e.g. ROID_TRIANGLE 	
			const ReferenceObjectID elemRoid = pElem->reference_object_id();
			
			try{
				//Do we need the type or use best?
				const QuadratureRule<dim>& rQuadRule = QuadratureRuleProvider<dim>::get(elemRoid, quadOrder);

			//	get reference element mapping by reference object id
				DimReferenceMapping<dim, WorldDim>& rMapping = ReferenceMappingProvider::
																get<dim, WorldDim>(elemRoid,vCorner);

			//	number of integration points
				const number* vWeight = rQuadRule.weights();
				const size_t numIP = rQuadRule.size();
				std::vector<MathVector<dim> > vLocalIP(numIP);
				std::vector<MathVector<dim> > vGlobalIP(numIP);
				
				//	Get Velocity Values in the Element
				//Stores integers with dof identifiers
				std::vector<DoFIndex> vInd;
				//vector of vectors sized #functions in function group
				std::vector<std::vector<number> > vvValue(vFctID.size());
				//We visit each vector inside the vvValue
				for(size_t fct = 0; fct < vvValue.size(); ++fct){
					spGridFct->dof_indices(pElem, vFctID[fct], vInd);
					vvValue[fct].resize(vInd.size());
					for(size_t sh = 0; sh < vInd.size(); ++sh){
						//DoFRef will give the value of GridFunction at index provided (dof)
						vvValue[fct][sh] = DoFRef(*spGridFct, vInd[sh]);
					}//end sh value for
				}//end fct value for

				
				for(size_t ip = 0; ip < numIP; ++ip)
					vLocalIP[ip]=rQuadRule.point(ip); 
				
				for(size_t ip = 0; ip < numIP; ++ip)
					rMapping.local_to_global(vGlobalIP[ip], vLocalIP[ip]);
				
				
				//Transformation matrices
				//std::vector<MathMatrix<dim, WorldDim> > vElemJT(numIP);
				vElemJT.resize(numIP);
				rMapping.jacobian_transposed(&(vElemJT[0]), &vLocalIP[0], numIP);
				
				//element integral value, reset to zero for each element
				number elemValue = 0;

				
				//loop integration points
				for(size_t ip = 0; ip < numIP; ++ip)
				{
					number ip_value=0.0;
					//	1. Interpolate values for each component inside a vector
					std::vector<number> vValuesAtIps(dim);
					for(int d1=0;d1<dim;d1++)
					{
						const LocalShapeFunctionSet<dim>& rTrialSpaceP1 =
														LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						std::vector<number> vShape;
						rTrialSpaceP1.shapes(vShape, vLocalIP[ip]);
					
						for(size_t sh = 0; sh < vvValue[d1].size(); ++sh){
							vValuesAtIps[d1] +=  vShape[sh]*vvValue[d1][sh];
						}
						
					}

					//  2. Interpolate gradient of the gridfunction
					//Stores the gradients (a vector size dim ) in a matrix with dim columns
					std::vector<MathVector<dim> > vvLocGradV[dim];
					std::vector<MathVector<dim> > vvGradV[dim];
					MathMatrix<dim, dim> JTInv;
					Inverse(JTInv, vElemJT[ip]);
					for(int d1 = 0; d1 < dim; ++d1){
						
						const LocalShapeFunctionSet<dim>& rTrialSpaceP = LocalFiniteElementProvider::get<dim>(elemRoid, vLFEID[d1]);
						// grads returns returns all gradients evaluated at a several points
						// This function returns the gradients of all Shape Functions at several element-local evaluation point in an array.
						//Parameters
						//	[out]	vvGrad	Vector of gradients
						//	[in]	vLocPos	Vector of Position on reference element

						rTrialSpaceP.grads(vvLocGradV[d1], vLocalIP[ip]);

						vvGradV[d1].resize(vvLocGradV[d1].size());
						for(size_t sh = 0; sh < vvGradV[d1].size(); ++sh){
							MatVecMult(vvGradV[d1][sh], JTInv, vvLocGradV[d1][sh]);
						}
								
					}

					MathMatrix<dim, dim> FunctionGradient;MatSet(FunctionGradient,0.0);
					for(int d1 = 0; d1 < dim; ++d1){
						for(int d2 = 0; d2 < dim; ++d2){
							FunctionGradient(d1, d2) = 0.0;
							for(size_t sh = 0; sh < vvValue[d1].size(); ++sh){							
								FunctionGradient(d1, d2) += vvValue[d1][sh] * vvGradV[d1][sh][d2];
							}
						}
					}
					
					//3. Perform the power functions for each term
					//We have to do: (u_p)^p + (Grad[u_p])^p
					for(int d1=0;d1<dim;d1++)
					{
						ip_value += std::pow(std::abs(vValuesAtIps[d1]), norm_order);
					}

					number mat_norm=0;
					for(int d1=0;d1<dim;d1++){
						for(int d2=0;d2<dim;d2++){
							mat_norm += std::pow(std::abs(FunctionGradient[d1][d2]),norm_order);
						}
					}
					ip_value += mat_norm;
					//5. Do the integration
					//	get quadrature weight
					const number weightIP = vWeight[ip];
					//const number weightIP = rQuadRule.weight(ip);
					//	get determinate of mapping
					const number det = SqrtGramDeterminant(vElemJT[ip]);
					
					elemValue += ip_value * weightIP * det ;					
				}//end ip for
				
				scalar_product += elemValue;	
				
			}UG_CATCH_THROW("VecProd failed.");//end try
		}//end element iteratior for
	}//end intSSGrp for
	//std::cout<<"L2VecProd: 5.- END OF LOOP FOR SUBSETS AND RETURN\n";
#ifdef UG_PARALLEL
	// sum over processes
	if(pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number local = scalar_product;
		com.allreduce(&local, &scalar_product, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
	}
#endif 
	//std::cout<<"L2VecProd: 6.- RESULT OF DOT PRODUCT: "<<scalar_product<<"\n";
	scalar_product=std::pow(scalar_product, 1.0/norm_order);
	return scalar_product;
	
}//end PNormToP

/// translate function names to IDs
template <typename TDomain, typename TAlgebra>
static size_t translateName2FctID(const GridFunction<TDomain, TAlgebra>& u,  const char* fctNames, std::vector<size_t> &vFct)
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




template <class TVector>
class ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero : public pcl::ICommunicationPolicy<IndexLayout>
{
	public:
	///	Default constructor
		ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero() : m_pVec(NULL)	{}

	///	Constructor setting the values
		ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero(TVector* pVec) 
	  	{
			set_vector(pVec);
		}

	///	sets the vector used in communication
		void set_vector(TVector* pVec) 
		{
			m_pVec = pVec;
			m_vProcessed.resize(m_pVec->size(), false);
		}

	/// clear processed flag
		void clear()
		{
			m_vProcessed.clear();
			m_vProcessed.resize(m_pVec->size(), false);
		}

	/// returns the buffer size
	/**
	 * This function returns the size of the buffer needed for the communication
	 * of passed interface. If the vector has fixed size entries this is just
	 * the number of interface entries times the size of the entry. In case
	 * of a variable size entry type a negative value is returned to indicate
	 * that no buffer size can be determined in advanced.
	 *
	 * \param[in]	interface	Interface that will communicate
	 */
		virtual int
		get_required_buffer_size(const Interface& interface)
		{
			if(block_traits<typename TVector::value_type>::is_static)
				return interface.size() * sizeof(typename TVector::value_type);
			else
				return -1;
		}

	///	writes the interface values into a buffer that will be sent
	/**
	 * This function collects all entries of the vector into a buffer that
	 * are part of the interface.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that will communicate
	 */
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{

			PROFILE_BEGIN_GROUP(ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero_collect, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			// get index
				const size_t index = interface.get_element(iter);

			// copy value
				Serialize(buff, v[index]);
				v[index] *= 0;
			}
			return true;
		}

	///	subtracts values of a buffer to the interface values
	/**
	 * This function subtracts the buffer values to the vector values.
	 *
	 * \param[out]		buff		Buffer
	 * \param[in]		interface	Interface that communicates
	 */
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			PROFILE_BEGIN_GROUP(ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero_extract, "algebra parallelization");
		//	check that vector has been set
			if(m_pVec == NULL) return false;

		//	rename for convenience
			TVector& v = *m_pVec;

		// entry
			typename TVector::value_type entry, diff;

		//	loop interface
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	get index
				const size_t index = interface.get_element(iter);

			//	copy vector
				Deserialize(buff, entry);

				if(v[index] == 0.0){
					if(entry != 0.0)
						v[index] = entry;
				} else {
					if(entry != 0.0)

						diff = v[index]; diff -= entry;
						if(diff != 0.0)
							UG_LOG("Should be the same entry, if nonzero.");
				}

			}
			return true;
		}

	private:
		TVector* m_pVec;
		std::vector<bool> m_vProcessed;
};


template <typename TVector>
void VecSetMasterToNonzeroValueIfAnyAndSlavesToZero(	TVector* pVec,
                            const IndexLayout& masterLayout,
                            const IndexLayout& slaveLayout,
                            pcl::InterfaceCommunicator<IndexLayout>* pCom = NULL)
{
	//	create a new communicator if required.
		pcl::InterfaceCommunicator<IndexLayout> tCom;
		if(!pCom)
			pCom = &tCom;
		pcl::InterfaceCommunicator<IndexLayout>& com = *pCom;

	//	create the required communication policies
		ComPol_VecSetMasterToNonzeroValueIfAnyAndSlavesToZero<TVector> cpVecSetNZifAny(pVec);

	//	sending: slaves, receiving: masters; masters subtract the value of only
	//	one slave on reception (according to the policy used)
		com.send_data(slaveLayout, cpVecSetNZifAny);
		com.receive_data(masterLayout, cpVecSetNZifAny);
		com.communicate();

}


/// set all vector values to zero for DoFs not connected to some (surface) subdomain 
template <typename TDomain, typename TAlgebra>
void SetZeroAwayFromSubset(	SmartPtr<GridFunction<TDomain, TAlgebra> > rhs, 
							const char* fctNames,
							const char* subsetNames)
{
    PROFILE_FUNC();

	if(!rhs->change_storage_type(PST_CONSISTENT))
		UG_THROW("SetZeroAwayFromSubset: Cannot change storage type");

	// dimension
	static const int dim = TDomain::dim;

	// translate name to ID
	std::vector<size_t> vFct(dim);
	translateName2FctID(*rhs, fctNames, vFct);

	// get subset to loop
	SubsetGroup ssGrp(rhs->domain()->subset_handler());
	ssGrp.add(TokenizeString(subsetNames));

	// get dof distribution
	SmartPtr<DoFDistribution> dd = rhs->dd();

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

		// check all neighbors
		std::vector<Vertex*> vNeighbors;
		
		CollectNeighbors(vNeighbors, *rhs->domain()->grid(), v, NHT_EDGE_NEIGHBORS);
		
		bool bConnected = false;
		for(size_t n = 0; n < vNeighbors.size(); ++n){
			int si = rhs->domain()->subset_handler()->get_subset_index(vNeighbors[n]);

			if(ssGrp.contains(si))
				bConnected = true;
		}
		if(bConnected)
			continue;

		//	index vector
		std::vector<DoFIndex> vMultInd;

		//	compute displacement (must loop all components)
		for(size_t i = 0; i < vFct.size(); ++i)
		{
			//	load indices associated with vertex
			if(dd->inner_dof_indices(v, vFct[i], vMultInd) == 1)
				DoFRef(*rhs, vMultInd[0]) = 0.0;
		}

	}

	
	VecSetMasterToNonzeroValueIfAnyAndSlavesToZero(&(*rhs), rhs->layouts()->master(), rhs->layouts()->slave(), &(rhs->layouts()->comm()) );
	
	rhs->set_storage_type(PST_ADDITIVE);

}


/// Modify domain (vertex positions) by displacement vector
template <typename TDomain, typename TAlgebra>
void TransformDomainByDisplacement(/*SmartPtr<TDomain> dom,*/
								   SmartPtr<GridFunction<TDomain, TAlgebra> > u,
								   const char* fctNames)
{
        PROFILE_FUNC();
	// dimension
	static const int dim = TDomain::dim;

	// translate name to ID
	std::vector<size_t> vFct(dim);
	translateName2FctID(*u, fctNames, vFct);


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
				"TransformDomainByDisplacement: No shift or shift with correct dimension!");
		// if (alpha == dim) std::cout << myVtxDisp << std::endl;


		// add displacement to vertex coordinates (on all levels)
		const MultiGrid& mg = *domain->grid();
		while(v){
			VecAppend(aaPos[v], myVtxDisp);
			//GridObject* parent = mg.get_parent(v);
			//if(parent)
			v = dynamic_cast<Vertex*>(mg.get_parent(v));
			// v=0x0, iff parent does not exist or not of vertex-type
		}

	}
}




/**
 * Adopted from InterpolateOnElements in interpolate.h
 * */

number CalculateTetrahedronVolume(const vector1& a, const vector1& b,
								  const vector1& c, const vector1& d)
{return 0.0;}

number CalculateTetrahedronVolume(const vector2& a, const vector2& b,
								  const vector2& c, const vector2& d)
{return 0.0;}

template <typename TGridFunction>
void CheckElementsForDisplacements(SmartPtr<TGridFunction> spGridFct, std::vector<size_t> &vFct)
{

	//	get reference element type
	typedef typename reference_element_traits<Tetrahedron>::reference_element_type
			ref_elem_type;
	const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

	//	dimension of reference element
	const int dim = ref_elem_type::dim;

	//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

	typename domain_type::grid_type &grid = *(spGridFct->domain()->grid());

	//	get iterators
	typename TGridFunction::template traits<Tetrahedron>::const_iterator iterEnd, iter;
	iterEnd = spGridFct->template end<Tetrahedron>(SurfaceView::ALL);
	iter = spGridFct->template begin<Tetrahedron>(SurfaceView::ALL);

	//	check if something to do:
	if(iter == iterEnd) return;


	//	id of shape functions used
    LFEID id = spGridFct->local_finite_element_id(vFct[0]);

	//	get trial space
    const LocalShapeFunctionSet<dim>& trialSpace =
    LocalFiniteElementProvider::get<dim>(roid, id);

	//	number of dofs on element
    const size_t nsh = trialSpace.num_sh();

    
    double volMin=1.0;
    double volMax=1.0;
    double elemDispMax=0.0;



	//	loop all elements
	for(;iter != iterEnd; ++iter)
	{
		double cornerDispMax=0.0;

		//	get element
		Tetrahedron* elem = *iter;

		std::vector<Vertex*> vElemVertices;
		CollectVertices(vElemVertices, grid, elem, true);

		//	determine CURRENT corner coordinates
		std::vector<position_type> vCornerOld;
		CollectCornerCoordinates(vCornerOld, *elem, *spGridFct->domain());


		// determine NEW corner coordinates
		std::vector<position_type> vCornerNew(nsh);
		UG_ASSERT(vCornerOld.size() == nsh, "Number of corners must match num shapes!")

		// 	loop all shapes
		for(size_t i = 0; i < nsh; ++i)
		{

			//  extract displacement
			position_type cornerDisp(0.0);

			//	index vector
			std::vector<DoFIndex> vMultInd;

			//	compute displacement (must loop all components)
			size_t alpha = 0;
			for(size_t a = 0; a < vFct.size(); ++a)
			{
				//	load indices associated with vertex
				if (spGridFct->dof_indices(vElemVertices[i], vFct[a], vMultInd)==1)
				{
					cornerDisp[alpha++] = DoFRef(*spGridFct, vMultInd[0]);
				}
			}
			UG_ASSERT(((alpha ==0) || (alpha==dim)),
					"CheckElementsForDisplacements: No shift or shift with correct dimension!");

			VecSubtract(vCornerNew[i], vCornerOld[i], cornerDisp);

			cornerDispMax = std::max(cornerDispMax,VecLength(cornerDisp));
		}

		// eval influence of displacement
		number vOld = CalculateTetrahedronVolume(vCornerOld[0], vCornerOld[1], vCornerOld[2], vCornerOld[3]);
		number vNew = CalculateTetrahedronVolume(vCornerNew[0], vCornerNew[1], vCornerNew[2], vCornerNew[3]);

        number ratio =vNew/vOld;
        volMin = std::min(ratio, volMin);
        volMax= std::max(ratio, volMax);



        elemDispMax = std::max(elemDispMax, cornerDispMax/pow(vOld, 1.0/3.0));


	}
#ifdef UG_PARALLEL
	ConstSmartPtr< AlgebraLayouts> layouts = spGridFct->dd()->layouts();
	volMin = layouts->proc_comm().allreduce(volMin, PCL_RO_MIN);
	volMax = layouts->proc_comm().allreduce(volMax, PCL_RO_MAX);
	elemDispMax = layouts->proc_comm().allreduce(elemDispMax, PCL_RO_MAX);
#endif
	UG_LOG ("Volumenaenderung: " << volMin << " --- "<< volMax<< std::endl);
	UG_LOG ("Max.Verschiebung: " << elemDispMax << std::endl);
}

/// Modify domain (vertex positions) by displacement vector
template <typename TDomain, typename TAlgebra>
void CheckDomainForDisplacements(SmartPtr<GridFunction<TDomain, TAlgebra> > u,
								  const char* fctNames)
{
	// translate names to IDs
	static const int dim = TDomain::dim;
	std::vector<size_t> vFct(dim);
	translateName2FctID(*u, fctNames, vFct);

	typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	CheckElementsForDisplacements<grid_function_type>(u, vFct);

}


template <typename TDomain, typename TAlgebra, typename TAlgebraSrc>
void AssignGridFunction(SmartPtr<GridFunction<TDomain, TAlgebra> > dst, const char* dstFctName,
		SmartPtr<GridFunction<TDomain, TAlgebraSrc> > src, const char* srcFctName, number scale=1.0, int si=-1)
{
	//	types
		typedef GridFunction<TDomain, TAlgebra> TFunction;
		typedef typename TFunction::domain_type::grid_type grid_type;
		typedef typename TFunction::element_type element_type;
//		const int dim = TDomain::dim;

	//	grid functions and function id
		GridFunction<TDomain, TAlgebraSrc> &usrc = (*src);
		const size_t srcid = usrc.fct_id_by_name(srcFctName);

		GridFunction<TDomain, TAlgebra> &udst = (*dst);
		const size_t dstid = udst.fct_id_by_name(dstFctName);


		// get dof distribution and index vectors
		SmartPtr<DoFDistribution> dstdist = udst.dd();
		SmartPtr<DoFDistribution> srcdist = usrc.dd();
		std::vector<DoFIndex> vSrcInd;
		std::vector<DoFIndex> vDstInd;


		typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		if (si<0)
		{
			// copy data for all vertices
			iter = dstdist->begin<Vertex>(SurfaceView::ALL); // SurfaceView::MG_ALL
			iterEnd = dstdist->end<Vertex>(SurfaceView::ALL);
		} else
		{
			// copy data only for subset si
			iter = dstdist->begin<Vertex>(si, SurfaceView::ALL); // SurfaceView::MG_ALL
			iterEnd = dstdist->end<Vertex>(si, SurfaceView::ALL);

		}


		// copy data
		for(;iter != iterEnd; ++iter)
		{
				Vertex* v = *iter;
				const size_t nfct = srcdist->inner_dof_indices(v, srcid, vSrcInd) ;
				if(dstdist->inner_dof_indices(v, dstid, vDstInd)==nfct)
				{

					for (size_t i=0; i<nfct; ++i)
					{
						number val = DoFRef(usrc, vSrcInd[i]);
						DoFRef(udst, vDstInd[i]) = scale*val;
					}
				}

		}

}


template <typename TDomain, typename TAlgebra, typename TAlgebraSrc>
void AddGridFunction(SmartPtr<GridFunction<TDomain, TAlgebra> > dst, const char* dstFctName,
		SmartPtr<GridFunction<TDomain, TAlgebraSrc> > src, const char* srcFctName, number scale=1.0, int si=-1)
{
	//	types
		typedef GridFunction<TDomain, TAlgebra> TFunction;
		typedef typename TFunction::domain_type::grid_type grid_type;
		typedef typename TFunction::element_type element_type;
//		const int dim = TDomain::dim;

	//	grid functions and function id
		GridFunction<TDomain, TAlgebraSrc> &usrc = (*src);
		const size_t srcid = usrc.fct_id_by_name(srcFctName);

		GridFunction<TDomain, TAlgebra> &udst = (*dst);
		const size_t dstid = udst.fct_id_by_name(dstFctName);


		// get dof distribution and index vectors
		SmartPtr<DoFDistribution> dstdist = udst.dd();
		SmartPtr<DoFDistribution> srcdist = usrc.dd();
		std::vector<DoFIndex> vSrcInd;
		std::vector<DoFIndex> vDstInd;


		typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;
		if (si<0)
		{
			// copy data for all vertices
			iter = dstdist->begin<Vertex>(SurfaceView::ALL); // SurfaceView::MG_ALL
			iterEnd = dstdist->end<Vertex>(SurfaceView::ALL);
		} else
		{
			// copy data only for subset si
			iter = dstdist->begin<Vertex>(si, SurfaceView::ALL); // SurfaceView::MG_ALL
			iterEnd = dstdist->end<Vertex>(si, SurfaceView::ALL);

		}


		// copy data
		for(;iter != iterEnd; ++iter)
		{
				Vertex* v = *iter;
				const size_t nfct = srcdist->inner_dof_indices(v, srcid, vSrcInd) ;
				if(dstdist->inner_dof_indices(v, dstid, vDstInd)==nfct)
				{

					for (size_t i=0; i<nfct; ++i)
					{
						number val = DoFRef(usrc, vSrcInd[i]);
						DoFRef(udst, vDstInd[i]) += scale*val;
					}
				}

		}

}


template <typename TDomain, typename TAlgebra>
void ClearGridFunction(SmartPtr<GridFunction<TDomain, TAlgebra> > dst, const char* dstFctNames, int si)
{
	//	types
		typedef GridFunction<TDomain, TAlgebra> TFunction;
		typedef typename TFunction::domain_type::grid_type grid_type;
		typedef typename TFunction::element_type element_type;
		const int dim = TDomain::dim;

		// ignore negative subsets
		if (si<0) return;

		//	grid functions and function ids
		GridFunction<TDomain, TAlgebra> &udst = (*dst);
		std::vector<size_t> vFct(dim);
		translateName2FctID(udst, dstFctNames, vFct);
		//const size_t dstid = udst.fct_id_by_name(dstFctNames);

		// get dof distribution and index vectors
		SmartPtr<DoFDistribution> dstdist = udst.dd();
		std::vector<DoFIndex> vDstInd;

		// clear data only for subset si
		typename DoFDistribution::traits<Vertex>::const_iterator iter, iterEnd;

		iter = dstdist->begin<Vertex>(si, SurfaceView::ALL); // SurfaceView::MG_ALL
		iterEnd = dstdist->end<Vertex>(si, SurfaceView::ALL);

		// clear data
		for(;iter != iterEnd; ++iter)
		{
			Vertex* v = *iter;

			for(size_t a = 0; a < vFct.size(); ++a)
			{
				const size_t nfct = dstdist->inner_dof_indices(v, vFct[a], vDstInd);

				for (size_t i=0; i<nfct; ++i)
				{ DoFRef(udst, vDstInd[i]) = 0.0; }
			}

		}

}


//! Assigns
template <typename TDomain, typename TAlgebra>
void AssignParallelPartitionOfUnity(SmartPtr<GridFunction<TDomain, TAlgebra> > dst, const char* dstFctName)
{

		GridFunction<TDomain, TAlgebra> &u = (*dst);
		u.set(1.0);

#ifdef UG_PARALLEL
		ConsistentToUnique(&u, u.layouts()->slave());
		u.set_storage_type(PST_ADDITIVE);
#endif
}


/// compute element area
/**
 * This is a test for the next function...
 */
template <typename TDomain, typename TAlgebra>
void AssignSurfElemArea(SmartPtr<GridFunction<TDomain, TAlgebra> > uSurf, const char* fctNameSurf, int maxLevel)
{

	//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	//typedef typename element_type::face_type face_type;
//	const int dim = TDomain::dim;

//	function id
	GridFunction<TDomain, TAlgebra> &u = (*uSurf);
	const size_t fctid = u.fct_id_by_name(fctNameSurf);

//	get multigrid
	if (u.domain()->grid().invalid()) return;
	//grid_type &grid = *(u.domain()->grid());

	// get dof distribution
	SmartPtr<DoFDistribution> dd = u.dd();

	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();

// collect adjacent volumes
	std::vector<Volume*> vAdjVolumes;

//	index vector
	std::vector<DoFIndex> vMultInd;

// evaluate area for all surface elements
	typename DoFDistribution::traits<Face>::const_iterator iter, iterEnd;
	iter = dd->begin<Face>(SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<Face>(SurfaceView::ALL);

	for(;iter != iterEnd; ++iter)
	{
		//	get face area
		Face* f = *iter;
		number area = FaceArea(f, aaPos);
		if(dd->inner_dof_indices(f, fctid, vMultInd) == 1)
			DoFRef(u, vMultInd[0]) = area;
	}

};


template <typename TAAPos>
void MyCalculateNormal(vector1& vNormOut, FaceVertices* face, TAAPos& aaPos)
{ return; }

template <typename TAAPos>
void MyCalculateNormal(vector2& vNormOut, FaceVertices* face, TAAPos& aaPos)
{ return; }

template <typename TAAPos>
void MyCalculateNormal(vector3& vNormOut, FaceVertices* face, TAAPos& aaPos)
{ CalculateNormal(vNormOut, face, aaPos); }

template <typename TAAPos>
void MyCalculateNormal(vector1& vNormOut, EdgeVertices* edge, TAAPos& aaPos)
{ return; }

template <typename TAAPos>
void MyCalculateNormal(vector2& vNormOut, EdgeVertices* edge, TAAPos& aaPos)
{ vNormOut=CalculateNormal(edge,aaPos); }

template <typename TAAPos>
void MyCalculateNormal(vector3& vNormOut, EdgeVertices* edge, TAAPos& aaPos)
{ return; }

template <typename TDomain, typename TAlgebra>
void AssignSurfElemNormal2d(GridFunction<TDomain, TAlgebra> &u, std::vector<size_t> &vFct, typename GridFunction<TDomain, TAlgebra>::domain_type::grid_type &grid, SmartPtr<DoFDistribution> dd, typename TDomain::position_accessor_type& aaPos, int surfid, number scale)
{
	const int dim = TDomain::dim;

	// collect adjacent vertices
	std::vector<Vertex*> vAdjVertices;


//	index vector
	std::vector<DoFIndex> vMultInd;

// evaluate area for all surface elements
	typename DoFDistribution::traits<Edge>::const_iterator iter, iterEnd;
	iter = dd->begin<Edge>(surfid, SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<Edge>(surfid, SurfaceView::ALL);

	MathVector<dim, number> normal;
	for(;iter != iterEnd; ++iter)
	{
		//	get face normal
		Edge* e = *iter;
		MyCalculateNormal(normal,e, aaPos);

		CollectVertices(vAdjVertices, grid, e);

		// assign to dofs
		size_t alpha=0;
		for(size_t i = 0; i < vFct.size(); ++i)
		{

			// a) assign for dofs on face
			if(dd->inner_dof_indices(e, vFct[i], vMultInd) == 1)
			{ DoFRef(u, vMultInd[0]) = normal[alpha]*scale; }


			// b) distribute for associated vertices
//			number vertCounter=0;
			for(std::vector<Vertex*>::iterator vert = vAdjVertices.begin(); vert != vAdjVertices.end(); ++vert)
			{

				// UG_LOG("Distributing" << dd->inner_dof_indices(*it, vFct[i], vMultInd));
				if(dd->inner_dof_indices(*vert, vFct[i], vMultInd) == 1)
				{ DoFRef(u, vMultInd[0]) += normal[alpha]*scale; }  // normalize
			}

			alpha++;
		}
	}

};

template <typename TDomain, typename TAlgebra>
void AssignSurfElemNormal3d(GridFunction<TDomain, TAlgebra> &u, std::vector<size_t> &vFct, typename GridFunction<TDomain, TAlgebra>::domain_type::grid_type &grid, SmartPtr<DoFDistribution> dd, typename TDomain::position_accessor_type& aaPos, int surfid, number scale)
{
	const int dim = TDomain::dim;

	// collect adjacent vertices
	std::vector<Vertex*> vAdjVertices;


//	index vector
	std::vector<DoFIndex> vMultInd;

	typename DoFDistribution::traits<Face>::const_iterator iter, iterEnd;
	iter = dd->begin<Face>(surfid, SurfaceView::ALL); // SurfaceView::MG_ALL
	iterEnd = dd->end<Face>(surfid, SurfaceView::ALL);

	MathVector<dim, number> normal;
	for(;iter != iterEnd; ++iter)
	{
		//	get face normal
		Face* f = *iter;
		MyCalculateNormal(normal,f, aaPos);

		CollectVertices(vAdjVertices, grid, f);

		// assign to dofs
		size_t alpha=0;
		for(size_t i = 0; i < vFct.size(); ++i)
		{

			// a) assign for dofs on face
			if(dd->inner_dof_indices(f, vFct[i], vMultInd) == 1)
			{ DoFRef(u, vMultInd[0]) = normal[alpha]*scale; }


			// b) distribute for associated vertices
			for(std::vector<Vertex*>::iterator vert = vAdjVertices.begin(); vert != vAdjVertices.end(); ++vert)
			{

				// UG_LOG("Distributing" << dd->inner_dof_indices(*it, vFct[i], vMultInd));
				if(dd->inner_dof_indices(*vert, vFct[i], vMultInd) == 1)
				{ DoFRef(u, vMultInd[0]) += normal[alpha]*scale;}  // normalize
			}
			alpha++;
		}
	}

};

template <typename TDomain, typename TAlgebra>
void AssignSurfElemNormal(SmartPtr<GridFunction<TDomain, TAlgebra> > uSurf, const char* fctNamesSurf, int surfid, number scale=1.0)
{

	//	types
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;
	const int dim = TDomain::dim;

//	function id
	GridFunction<TDomain, TAlgebra> &u = (*uSurf);
	//const size_t fctid = u.fct_id_by_name(fctNameSurf);

	// translate name to ID
	std::vector<size_t> vFct(dim);
	translateName2FctID(u, fctNamesSurf, vFct);


//	get multigrid
	if (u.domain()->grid().invalid()) return;
	grid_type &grid = *(u.domain()->grid());

	// get dof distribution
	SmartPtr<DoFDistribution> dd = u.dd();

	//	get position accessor
	typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();

	if(dim==2){
		AssignSurfElemNormal2d(u, vFct, grid, dd, aaPos, surfid, scale);
	}
	if(dim==3){
		AssignSurfElemNormal3d(u, vFct, grid, dd, aaPos, surfid, scale);
	}

#ifdef UG_PARALLEL
	//UG_LOG("ADDITIVE")
    // changing back from additive to consistent
   // uSurf->set_storage_type(PST_ADDITIVE);
#endif


};

template <typename TDomain, typename TAlgebra>
void AssignSurfElemNormalFV(SmartPtr<GridFunction<TDomain, TAlgebra> > uPtr, const char* fctNames,
                                   const char* BndSubset, const char* InnerSubset, std::vector<number> BndOrientation)
{

//	types
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

	GridFunction<TDomain, TAlgebra> &u = (*uPtr);
//	read subsets
	SubsetGroup innerSSGrp(u.domain()->subset_handler());
	if(InnerSubset != NULL)
		innerSSGrp.add(TokenizeString(InnerSubset));
	else // add all if no subset specified
		innerSSGrp.add_all();

//	read bnd subsets
	SubsetGroup bndSSGrp(u.domain()->subset_handler());
	if(BndSubset != NULL){
		bndSSGrp.add(TokenizeString(BndSubset));
	}
	else{
		UG_THROW("ShapeOptim::AssignSurfElemNormalFV: No boundary subsets specified. Aborting.");
	}

	//	loop subsets
	for(size_t i = 0; i < innerSSGrp.size(); ++i)
	{
		//	get subset index
		const int si = innerSSGrp[i];

		if (innerSSGrp.dim(i) != TGridFunction::dim)
			UG_THROW("ShapeOptim::AssignSurfElemNormalFV: Element dimension does not match world dimension!");

		//	create integration kernel
		static const int dim = TGridFunction::dim;

		// translate name to ID
		std::vector<size_t> vFct(dim);
		translateName2FctID(u, fctNames, vFct);

		//	integrate elements of subset
		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

		//	get iterators for all elems on subset
		typedef typename TGridFunction::template dim_traits<dim>::const_iterator const_iterator;
		const_iterator iter = u.template begin<grid_base_object>(si);
		const_iterator iterEnd = u.template end<grid_base_object>(si);

		//	create a FV1 Geometry
		DimFV1Geometry<dim> geo;

		//	specify, which subsets are boundary
		for(size_t s = 0; s < bndSSGrp.size(); ++s)
		{
			//	get subset index
			const int bndSubset = bndSSGrp[s];

			//	request this subset index as boundary subset. This will force the
			//	creation of boundary subsets when calling geo.update
			geo.add_boundary_subset(bndSubset);
		}

		//	vector of corner coordinates of element corners (to be filled for each elem)
		std::vector<MathVector<dim> > vCorner;

		//	loop elements of subset
		for( ; iter != iterEnd; ++iter)
		{
			//	get element
			grid_base_object* elem = *iter;

			//	get all corner coordinates
			CollectCornerCoordinates(vCorner, *elem, u.domain()->position_accessor(), true);

			//	compute bf at bip for element
			try{
				geo.update(elem, &vCorner[0], u.domain()->subset_handler().get());
			}
			UG_CATCH_THROW("ShapeOptim::AssignSurfElemNormalFV: "
							"Cannot update Finite Volume Geometry.");

			//	specify, which subsets are boundary
			for(size_t s = 0; s < bndSSGrp.size(); ++s)
			{
				//	get subset index
				const int bndSubset = bndSSGrp[s];
				const float bndOrient = BndOrientation[s];

				//	get all bf of this subset
				typedef typename DimFV1Geometry<dim>::BF BF;
				const std::vector<BF>& vBF = geo.bf(bndSubset);

				//	loop boundary faces
				for(size_t b = 0; b < vBF.size(); ++b)
				{
					//	get bf
					const BF& bf = vBF[b];

					//	get normal on bf
					const MathVector<dim>& normal = bf.normal();
					MathVector<dim> n;
					VecNormalize(n,normal);

					// assign to dofs
					for(size_t fi = 0; fi < vFct.size(); ++fi)
					{
						//	skip if function is not defined in subset
						if(!u.is_def_in_subset(vFct[fi], si)) continue;

						//	get fct multi-indices of element
						std::vector<DoFIndex> ind;
						u.dof_indices(elem, vFct[fi], ind);

						//	check multi indices
						UG_ASSERT(ind.size() == bf.num_sh(),
								"ShapeOptim::AssignSurfElemNormalFV: Wrong number of"
								" multi indices, ind: "<<ind.size() << ", bf.num_sh: "
								<< bf.num_sh());

						// 	assign normal to vertices (= degrees of freedom)
						for(size_t sh = 0; sh < bf.num_sh(); ++sh)
						{
							DoFRef(u, ind[sh]) += bndOrient*n[fi]*bf.shape(sh);
						}
					}
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra, typename TFunctor>
void AssignSurfFunction(SmartPtr<GridFunction<TDomain, TAlgebra> > uSurf, const char* fctNameSurf,
				 /* GridFunction<TDomain, TAlgebra>& uVol, const char* fctNameVol,*/
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userGradP,
				int si)
{
  PROFILE_FUNC();
//	types
	static const size_t dim = TDomain::dim;
	typedef GridFunction<TDomain, TAlgebra> TFunction;
	typedef typename TFunction::domain_type::grid_type grid_type;
	typedef typename TFunction::element_type element_type;

	//if(!importFlowY.data_given())
	//{ std::cout << "Missing import..." << std::endl; return; }

	// translate name to IDs
	GridFunction<TDomain, TAlgebra> &u = (*uSurf);
	std::vector<size_t> vFct(dim);
	if (translateName2FctID(u, fctNameSurf, vFct) == 0) return;

//	get multigrid
	if (u.domain()->grid().invalid()) return;

	// get dof distribution
	ConstSmartPtr<DoFDistribution> dd = u.dd();

#ifdef UG_PARALLEL
    // we will write values in an additive fashion (i.e. locally on each proc)
    uSurf->set_storage_type(PST_ADDITIVE);
#endif

    // evaluate
	//typename DoFDistribution::traits<Face>::const_iterator iter, iterEnd;
	//iter = u.template begin<Face>(SurfaceView::ALL); // SurfaceView::MG_ALL
	//iterEnd = u.template end<Face>(SurfaceView::ALL);

	typedef typename IteratorProvider<TFunction>::template traits<Face>::const_iterator const_iterator;
	const_iterator iter = IteratorProvider<TFunction>::template begin<Face>(u, si);
	const_iterator iterEnd = IteratorProvider<TFunction>::template end<Face>(u, si);


	int nFace = 0;
//	int nFaceUsed = 0;
	for(;iter != iterEnd; ++iter)
	{
		//	get face
		Face* face = *iter;

		// read data from functor
		typename TFunctor::return_type fvalue;
		TFunctor::apply(face, u, userFlowY, userGradP, fvalue);

		// assign data to vector
		size_t alpha=0;
		std::vector<DoFIndex> vMultInd;

		for(size_t i = 0; i < vFct.size(); ++i)
		{
			size_t ncmpi = dd->inner_dof_indices(face, vFct[i], vMultInd);
            for(size_t j=0; j<ncmpi; j++)
            { DoFRef(u, vMultInd[j]) = fvalue[alpha++]; }
		}

		UG_ASSERT(vFct.size()==TFunctor::return_size, "Sizes do not match (" << vFct.size() << "!=" << TFunctor::return_size << ")");
		UG_ASSERT((alpha==TFunctor::return_size || alpha == 0), "Sizes do not match (" << alpha << "!=" << TFunctor::return_size << ")");
		 nFace++;
	}
	UG_LOG("Iterated over "<< nFace << " faces,")


#ifdef UG_PARALLEL
    // changing back from additive to consistent 
    uSurf->change_storage_type(PST_CONSISTENT);
#endif
};




template <int dim>
struct SurfaceEnergyAvgFunctor{

	typedef number return_type;
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
inline	static double apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP)
	{ return 0.0; }

};

template <>
struct SurfaceEnergyAvgFunctor<3>{

	typedef number return_type;
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
	inline static double apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP)
		{
			typedef GridFunction<TDomain, TAlgebra> TFunction;

			//	get position accessor
			typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();
			typename TDomain::grid_type &grid = *(u.domain()->grid());

			// collect adjacent volumes
				std::vector<Volume*> vAdjVolumes;

			CollectVolumes(vAdjVolumes, grid, f);
			int nVolumes=0;
			number delta = 0.0;
			MathVector<3> fnormal3; // do this in 3D
			CalculateNormal(fnormal3, f, aaPos);



			for(size_t j = 0; j < vAdjVolumes.size(); ++j)
			{

				static const size_t dim = TDomain::dim;

				// Evaluate element
				Volume *elem = vAdjVolumes[j];
				const size_t numCo = elem->num_vertices();
				std::vector<MathVector<dim> > vCorner(numCo);
				CollectCornerCoordinates(vCorner, *elem, aaPos, true);

				// Evaluate gradient at element mid-point
				MathVector<dim> importVecY, importVecP;
				MathVector<dim> globIP;
				{
					// set integration point
					MathVector<dim> localIP(0.25);
					AveragePositions(globIP, &vCorner[0], numCo);

					//	create storage
					LocalIndices ind;
					LocalVector locU;

					// 	get global indices
					u.indices(elem, ind);

					// 	adapt local algebra
					locU.resize(ind);

					// 	read local values of u
					GetLocalVector(locU, u);
					int theSI = u.domain()->subset_handler()->get_subset_index(elem);


					//	compute data
					try{
						(*userFlowY)(importVecY, globIP, 0.0, theSI, elem,
								&vCorner[0], localIP, &locU);
						(*userFlowP)(importVecP, globIP, 0.0, theSI, elem,
								&vCorner[0], localIP, &locU);
					}
					UG_CATCH_THROW("EnergyAvg: Cannot evaluate data.");
				}


				delta += VecDot(importVecY, importVecP);
				nVolumes++;
			}

			return delta;
		}


};



template <int dim>
struct SurfaceEnergyJumpFunctor{

	typedef number return_type;
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
	static bool apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
						SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
						SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP,
						return_type &value)
	{ value = 0.0; return true;}

};

template <>
struct SurfaceEnergyJumpFunctor<3>{

	typedef number return_type;
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
	static double apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP)
		{
			typedef GridFunction<TDomain, TAlgebra> TFunction;

			//	get position accessor
			typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();
			typename TDomain::grid_type &grid = *(u.domain()->grid());

			// collect adjacent volumes
				std::vector<Volume*> vAdjVolumes;

			CollectVolumes(vAdjVolumes, grid, f);
			int nVolumes=0;
			number delta = 0.0;
			MathVector<3> fnormal3; // do this in 3D
			CalculateNormal(fnormal3, f, aaPos);


			for(size_t j = 0; j < vAdjVolumes.size(); ++j)
			{

				static const size_t dim = TDomain::dim;

				// Evaluate element
				Volume *elem = vAdjVolumes[j];
				const size_t numCo = elem->num_vertices();
				std::vector<MathVector<dim> > vCorner(numCo);
				CollectCornerCoordinates(vCorner, *elem, aaPos, true);

				// Evaluate gradient at element mid-point
				MathVector<dim> importVecY, importVecP;
				MathVector<dim> globIP;
				{
					// set integration point
					MathVector<dim> localIP(0.5);
					AveragePositions(globIP, &vCorner[0], numCo);

					//	create storage
					LocalIndices ind;
					LocalVector locU;

					// 	get global indices
					u.indices(elem, ind);

					// 	adapt local algebra
					locU.resize(ind);

					// 	read local values of u
					GetLocalVector(locU, u);
					int theSI = u.domain()->subset_handler()->get_subset_index(elem);


					//	compute data
					try{
						(*userFlowY)(importVecY, globIP, 0.0, theSI, elem,
								&vCorner[0], localIP, &locU);
						(*userFlowP)(importVecP, globIP, 0.0, theSI, elem,
								&vCorner[0], localIP, &locU);
					}
					UG_CATCH_THROW("EnergyJump: Cannot evaluate data.");
				}


				// Evaluate element orientation wrt face
				double angle;
				{
					MathVector<dim> fdelta(globIP);
					CollectCornerCoordinates(vCorner, f[0], aaPos, true);
					VecScaleAppend(fdelta, -1.0, vCorner[0]);

					MathVector<3> fdelta3;
					VecCopy(fdelta3, fdelta, 0.0);
					angle=VecDot(fnormal3, fdelta3);
				}


				delta += VecDot(importVecY, importVecP) * ((angle>0) - (angle<0));
				nVolumes++;
			}

			return delta;
		}

};









template <int dim>
struct SurfaceShapeGradValFunctor{

	typedef double return_type[1];
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
	static bool apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP,
					return_type &value)
	{ value[0] = 0.0; return true;}

};

template <>
struct SurfaceShapeGradValFunctor<3>{

	typedef double return_type[1];
	static size_t const return_size=1;

	template <typename TDomain, typename TAlgebra>
	  inline static bool apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userGradP,
					return_type &value)
			{
			  // PROFILE_FUNC();
				typedef GridFunction<TDomain, TAlgebra> TFunction;

				// get position accessor
				typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();
				typename TDomain::grid_type &grid = *(u.domain()->grid());

				// collect adjacent volumes
				std::vector<Volume*> vAdjVolumes;

				CollectVolumes(vAdjVolumes, grid, f);
				int nVolumes=0;
				number delta = 0.0;
				MathVector<3> fnormal3; // do this in 3D
				CalculateNormal(fnormal3, f, aaPos);

				for(size_t j = 0; j < vAdjVolumes.size(); ++j)
				{
					static const size_t dim = TDomain::dim;

					// Evaluate element
					Volume *elem = vAdjVolumes[j];
					const size_t numCo = elem->num_vertices();
					std::vector<MathVector<dim> > vCorner(numCo);
					CollectCornerCoordinates(vCorner, elem, aaPos, true);

					// Evaluate gradient at element mid-point
					MathVector<dim> importFlowY, importGradP;
					MathVector<dim> globIP;
					{
						// set integration point
						MathVector<dim> localIP(0.5);
						AveragePositions(globIP, &vCorner[0], numCo);

						//	create storage
						LocalIndices ind;
						LocalVector locU;

						// 	get global indices
						u.indices(elem, ind);

						// 	adapt local algebra
						locU.resize(ind);

						// 	read local values of u
						GetLocalVector(locU, u);
						int theSI = u.domain()->subset_handler()->get_subset_index(elem);


						//	compute data
						//		try{
							(*userFlowY)(importFlowY, globIP, 0.0, theSI, elem,
									&vCorner[0], localIP, &locU);
							(*userGradP)(importGradP, globIP, 0.0, theSI, elem,
									&vCorner[0], localIP, &locU);
							//}
							//UG_CATCH_THROW("ShapeGradValue: Cannot evaluate data.");
					}


					delta += VecDot(importFlowY, importGradP) - 2.0* VecDot(importFlowY, fnormal3)*VecDot(importGradP, fnormal3);
					nVolumes++;
				}

                UG_ASSERT(nVolumes==2,"ERROR: (Interior?) face does not have exactly 2 associated volumes!")

                value[0] = delta;
				return true;
			}
};








template <int dim>
struct SurfaceShapeGradFunctor{

	typedef MathVector<dim> return_type;
	static size_t const return_size=dim;

	template <typename TDomain, typename TAlgebra>
	static bool apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowP,
					return_type &value)
	{ value = 0.0; return true;}

};


//! Evaluates \vec n* [[ -2*(k \frac{\partial y}){\partial n}\frac{\partial p}{\partial n}  + (k \nabla y) \nabla p]]
template <>
struct SurfaceShapeGradFunctor<3>
{

	typedef MathVector<3> return_type;
	static size_t const return_size=3;

	template <typename TDomain, typename TAlgebra>
	static bool apply(Face *f, GridFunction<TDomain, TAlgebra> &u,
					  SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				      SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userGradP,
				      return_type &value)
	{ //PROFILE_FUNC();
		typedef GridFunction<TDomain, TAlgebra> TFunction;

		//	get position accessor
		typename TDomain::position_accessor_type& aaPos = u.domain()->position_accessor();
		typename TDomain::grid_type &grid = *(u.domain()->grid());

		// collect adjacent volumes
		std::vector<Volume*> vAdjVolumes;

		CollectVolumes(vAdjVolumes, grid, f);
		int nVolumes=0;
		number delta = 0.0;
		MathVector<3> fnormal3; // do this in 3D
		CalculateNormal(fnormal3, f, aaPos);

        double localCoeff;

		
        for(size_t j = 0; j < vAdjVolumes.size(); ++j)
		{

			static const size_t dim = TDomain::dim;

			// Evaluate element
			Volume *elem = vAdjVolumes[j];
			const size_t numCo = elem->num_vertices();
			std::vector<MathVector<dim> > vCorner(numCo);
			CollectCornerCoordinates(vCorner, *elem, aaPos, true);

			// Evaluate gradient at element mid-point
			MathVector<dim> importFlowY, importGradP;
			MathVector<dim> globIP;

		    // set integration point
		    MathVector<dim> localIP(0.25); // CHECK!!!
		    AveragePositions(globIP, &vCorner[0], numCo);

		    //	create storage
		    LocalIndices ind;
		    LocalVector locU;

		    // 	get global indices
		    u.indices(elem, ind);

		    // 	adapt local algebra
		    locU.resize(ind);

		    // 	read local values of u
		    GetLocalVector(locU, u);

            int isInner = 0;
		    const int   subsetInd  = u.domain()->subset_handler()->get_subset_index(elem);
            const char* subsetName = u.domain()->subset_handler()->get_subset_name(subsetInd);
            if( strcmp(subsetName, "V_INNER") )
            {
                isInner = 1;
            }
            else if( strcmp(subsetName, "V_OUTER") )
            {
                isInner = 0;
            }
            else
            {
                std::cout << "Subsets are broken" << std::endl;
                exit(1);
            }

		    //	compute data
//		    try{
			    (*userFlowY)(importFlowY, globIP, 0.0, subsetInd, elem,
					    &vCorner[0], localIP, &locU);
			    (*userGradP)(importGradP, globIP, 0.0, subsetInd, elem,
					    &vCorner[0], localIP, &locU);
//		    }
//		    UG_CATCH_THROW("ShapeGradValue: Cannot evaluate data.");

		    // Evaluate element orientation wrt face
		    double angle;
		    {
                CollectCornerCoordinates(vCorner, *elem, aaPos, true);
			    MathVector<dim> elemBaryCenter(globIP);
                AveragePositions(elemBaryCenter, &vCorner[0], elem->num_vertices());

			    CollectCornerCoordinates(vCorner, f[0], aaPos, true);
                MathVector<dim> faceBaryCenter(globIP);
                AveragePositions(faceBaryCenter, &vCorner[0], f->num_vertices());

			    MathVector<3> fdelta(elemBaryCenter);
			    VecScaleAppend(fdelta, -1.0, faceBaryCenter);

			    angle=VecDot(fnormal3, fdelta);

//                if(vAdjVolumes.size() != 2) std::cout << elemBaryCenter << std::endl;
		    }
            if (  isInner && angle > 0 ) VecScale(fnormal3, fnormal3, -1.0);
            if ( !isInner && angle < 0 ) VecScale(fnormal3, fnormal3, -1.0);

            // TODO: resolve this! Coefficients should be provided as arguments
            if(isInner) localCoeff = 1e-3;
            else localCoeff = 1.0;

			delta += (1.0 - 2.0*isInner)*localCoeff* ( VecDot(importFlowY, importGradP) - 2.0* VecDot(importFlowY, fnormal3)*VecDot(importGradP, fnormal3) );
			nVolumes++;
		} // AdjacentElements
		VecScale(value, fnormal3, delta);
		return true;
	}

};


template <typename TDomain, typename TAlgebra>
void AssignSurfaceShapeGradientValue(SmartPtr<GridFunction<TDomain, TAlgebra> > uSurf, const char* fctNameSurf,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userGradP,
				int si)
{
	AssignSurfFunction<TDomain, TAlgebra, SurfaceShapeGradValFunctor<TDomain::dim> >
	(uSurf, fctNameSurf, userFlowY, userGradP,si);
}


template <typename TDomain, typename TAlgebra>
void AssignSurfaceShapeGradient(SmartPtr<GridFunction<TDomain, TAlgebra> > uSurf, const char* fctNameSurf,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userFlowY,
				SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > userGradP,
				int si)
{
	AssignSurfFunction<TDomain, TAlgebra, SurfaceShapeGradFunctor<TDomain::dim> >
	(uSurf, fctNameSurf, userFlowY, userGradP,si);
}



} // namespace ShapeOptim
} // namespace ug
#endif
