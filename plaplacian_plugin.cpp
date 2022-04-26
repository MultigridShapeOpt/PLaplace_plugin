#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_algebra/operator/debug_writer.h"


#include "NavierStokes.h"
#include "NavierStokesAdjoint.h"
#include "PLaplaceDerivative.h"
#include "PLaplaceDerivativeRHS.h"
#include "VolumeConstraintSecondDerivative.h"
#include "XBarycenterConstraintSecondDerivative.h"
#include "plaplacian_functions.h"
#include "LargeProblemRHS.h"
#include "JPrime.h"
#include "TimeDependentNavierStokes.h"
#include "AdvancedGridFunction.h"

#include "restart_functions.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace PLaplacian{

struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();
//	TimeDependent NavierStokes Class with debugging functions, Non-linear
	{
		typedef TimeDependentNavierStokes<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("TimeDependentNavierStokes").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_nonlinear", &T::set_nonlinear, "", "Non-linear problem")
			.add_method("set_picard", &T::set_picard,"","Picard Iteration")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_stabilization_type", &T::set_stabilization_type)
			.add_method("set_time_delta", &T::set_time_delta)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"TimeDependentNavierStokes",tag);
	}

//	NavierStokes Class with debugging functions, Non-linear
	{
		typedef NavierStokes<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("IncompressibleNavierStokes").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_nonlinear", &T::set_nonlinear, "", "Non-linear problem")
			.add_method("set_picard", &T::set_picard,"","Picard Iteration")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_stabilization_type", &T::set_stabilization_type)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"IncompressibleNavierStokes",tag);
	}

//	AdjointNavierStokes Class with debugging functions, Non-linear
	{
		typedef NavierStokesAdjoint<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("NavierStokesAdjoint").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			//.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			//.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	//#ifdef UG_FOR_LUA
			//.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			//.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	//#endif
			.add_method("set_kinematic_viscosity", &T::set_kinematic_viscosity)
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_stabilization", &T::set_stabilization)
			.add_method("set_stabilization_type", &T::set_stabilization_type)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"NavierStokesAdjoint",tag);
	}
//	PLaplaceDerivative Class, linear, rhs used as defect of nonlinear problem, but checked externally
	{
		typedef PLaplaceDerivative<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("PLaplaceDerivative").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_lambda_vol", &T::set_lambda_vol)
			.add_method("set_lambda_barycenter", &T::set_lambda_barycenter)
			.add_method("set_p", &T::set_p)
			.add_method("set_control", &T::set_control)
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_epsilon", &T::set_epsilon)
			.add_method("set_elliptic_epsilon", &T::set_elliptic_epsilon)
			.add_method("set_step_length", &T::set_step_length)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"PLaplaceDerivative",tag);
	}
//	PLaplaceDerivativeRHS Class, linear, rhs used as defect of nonlinear problem, but checked externally
	{
		typedef PLaplaceDerivativeRHS<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("PLaplaceDerivativeRHS").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_pressure", static_cast<void (T::*)(number)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(number)>(&T::set_adjoint_pressure),"","AdjointPressureScalar")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_lambda_vol", &T::set_lambda_vol)
			.add_method("set_lambda_barycenter", &T::set_lambda_barycenter)
			.add_method("set_p", &T::set_p)
			.add_method("set_control", &T::set_control)
			.add_method("set_epsilon", &T::set_epsilon)
			.add_method("set_step_length", &T::set_step_length)
			.add_method("set_flag_sensitivity", &T::set_flag_sensitivity)
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"PLaplaceDerivativeRHS",tag);
	}
//	VolumeConstraintSecondDerivative 
	{
		typedef VolumeConstraintSecondDerivative<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("VolumeConstraintSecondDerivative").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_multiplier", &T::set_multiplier)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"VolumeConstraintSecondDerivative",tag);
	}
//	XBarycenterConstraintSecondDerivative 
	{
		typedef XBarycenterConstraintSecondDerivative<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("XBarycenterConstraintSecondDerivative").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_index", &T::set_index)
			.add_method("set_multiplier", &T::set_multiplier)
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"XBarycenterConstraintSecondDerivative",tag);
	}

//	RHS equal to Lu-B.delta_lambda	
	{
		typedef LargeProblemRHS<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("LargeProblemRHS").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d1", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d1),"","DeformationVectord1")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d2", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d2),"","DeformationVectord2")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_vector_d3", static_cast<void (T::*)(number)>(&T::set_deformation_vector_d3),"","DeformationVectord3")
			.add_method("set_deformation_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d1", static_cast<void (T::*)(number)>(&T::set_deformation_d1),"","DeformationValued1")
			.add_method("set_deformation_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d2", static_cast<void (T::*)(number)>(&T::set_deformation_d2),"","DeformationValued2")
			.add_method("set_deformation_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_deformation_d3", static_cast<void (T::*)(number)>(&T::set_deformation_d3),"","DeformationValued3")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")

			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_pressure", static_cast<void (T::*)(number)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(number)>(&T::set_adjoint_pressure),"","AdjointPressureScalar")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_lambda_vol", &T::set_lambda_vol)
			.add_method("set_lambda_barycenter", &T::set_lambda_barycenter)
			.add_method("set_p", &T::set_p)
			.add_method("set_epsilon", &T::set_epsilon)
			.add_method("set_control", &T::set_control)
			.add_method("set_step_length", &T::set_step_length)
			.add_method("set_multiplier_vol", &T::set_multiplier_vol)
			.add_method("set_multiplier_bx", &T::set_multiplier_bx)
			.add_method("set_multiplier_by", &T::set_multiplier_by)
			.add_method("set_multiplier_bz", &T::set_multiplier_bz)
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"LargeProblemRHS",tag);
	}

	//	Provides JPrime as a rhs
	{
		typedef JPrime<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("JPrime").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Functions#Subset(s)")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(number)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#ifdef UG_FOR_LUA
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(const char*)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
			.add_method("set_kinematic_viscosity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_kinematic_viscosity),"","KinematicViscosity")
	#endif
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d1),"","VelocityGradientd1")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")
			.add_method("set_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d2),"","VelocityGradientd2")

			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_velocity_vector_d3),"","VelocityGradientd3")
			.add_method("set_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d1", static_cast<void (T::*)(number)>(&T::set_velocity_d1),"","VelocityValued1")
			.add_method("set_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d2", static_cast<void (T::*)(number)>(&T::set_velocity_d2),"","VelocityValued2")
			.add_method("set_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_velocity_d3", static_cast<void (T::*)(number)>(&T::set_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d1),"","AdjointVelocityGradientd1")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d2),"","AdjointVelocityGradientd2")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>,dim>>)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_adjoint_velocity_vector_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_vector_d3),"","AdjointVelocityGradientd3")
			.add_method("set_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_pressure", static_cast<void (T::*)(number)>(&T::set_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_pressure),"","PressureScalar")
			.add_method("set_adjoint_pressure", static_cast<void (T::*)(number)>(&T::set_adjoint_pressure),"","AdjointPressureScalar")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d1", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d1),"","VelocityValued1")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d2", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d2),"","VelocityValued2")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")
			.add_method("set_adjoint_velocity_d3", static_cast<void (T::*)(number)>(&T::set_adjoint_velocity_d3),"","VelocityValued3")	
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name,"JPrime",tag);
	}
	
 
}

template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	
}

template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> approximation_space_type;
	typedef GridFunction<TDomain, TAlgebra> grid_func_type;
	//typedef AdvancedGridFunction<TDomain, TAlgebra> TDerived;
	//	AdvancedGridFunction
	{
		typedef AdvancedGridFunction<TDomain, TAlgebra> TDerived;
		string name = string("AdvancedGridFunction").append(suffix);
		reg.add_class_<TDerived, grid_func_type>(name, grp)
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>)>("ApproximationSpace")
			.template add_constructor<void (*)(SmartPtr<approximation_space_type>, int)>("ApproximationSpace#Level")
			.add_method("change_storage_type_to_consistent", &TDerived::change_storage_type_to_consistent)
			.add_method("change_storage_type_to_additive", &TDerived::change_storage_type_to_additive)
			.add_method("change_storage_type_to_unique", &TDerived::change_storage_type_to_unique)
			.add_method("has_storage_type_consistent", &TDerived::has_storage_type_consistent)
			.add_method("has_storage_type_additive", &TDerived::has_storage_type_additive)
			.add_method("has_storage_type_unique", &TDerived::has_storage_type_unique)
			.add_method("enforce_consistency", &TDerived::enforce_consistency)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "AdvancedGridFunction", tag);
	}

	reg.add_function("TransformDomainByDisplacement", static_cast<void (*) (SmartPtr<grid_func_type>, const char*) >
					(&ug::PLaplacian::TransformDomainByDisplacement<TDomain, TAlgebra>),
					grp,"","GridFunction#FunctionNames"
					);

	//Restart functions
	{
		reg.add_function("SaveNodalPositions2GridFunction", &SaveNodalPositions2GridFunction<TDomain,TAlgebra>, grp);
	}
	//Restart functions
	{
		reg.add_function("SetDomainCoordinatesFromGF", &SetDomainCoordinatesFromGF<TDomain,TAlgebra>, grp);
	}
	{
		reg.add_function("SetZeroAwayFromSubset",
							static_cast<void (*) (SmartPtr<grid_func_type>, const char*,  const char*) >
							(&ug::PLaplacian::SetZeroAwayFromSubset<TDomain, TAlgebra>),
						    grp, "", "GridFunction#FunctionNames#SubsetNames");
	}
	//P norm to the p
	{
		reg.add_function("PNormToP", &PNormToP<grid_func_type>, grp);
	}
}

}; // end Functionality

/**
 * Class exporting the functionality of the plugin restricted to 2 and 3 spatial
 * dimensions. All functionality that is to be used in scripts or visualization
 * only in 2d and 3d must be registered here.
 */

} 


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_PLaplacian(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/ElemDisc");
	typedef PLaplacian::Functionality Functionality;

	try{
		
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
