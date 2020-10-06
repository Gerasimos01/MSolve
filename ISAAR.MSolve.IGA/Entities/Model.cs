using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Transfer;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Transfer;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Elements.Boundary;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

[assembly: InternalsVisibleTo("ISAAR.MSolve.IGA.Tests")]

namespace ISAAR.MSolve.IGA.Entities
{
    /// <summary>
	/// Model that contains all data needed for Isogeometric analysis.
	/// </summary>
	public class Model : IModel
	{
		private IGlobalFreeDofOrdering _globalDofOrdering;


		private IList<PenaltyDofPair> PenaltyDofPairs { get; } = new List<PenaltyDofPair>();
		
		/// <summary>
		/// <see cref="Table{TRow,TColumn,TValue}"/> that contains the constrained degree of freedom and the value of their constrains.
		/// </summary>
		public Table<INode, IDofType, double> Constraints { get; private set; } =
			new Table<INode, IDofType, double>();

		/// <summary>
		/// Return an <see cref="IEnumerable{ControlPoint}"/> with the Control Points of the <see cref="Model"/>.
		/// </summary>
		public IEnumerable<ControlPoint> ControlPoints => ControlPointsDictionary.Values;

		/// <summary>
		/// Dictionary containing the Control Points of the model.
		/// </summary>
		public Dictionary<int, ControlPoint> ControlPointsDictionary { get; } = new Dictionary<int, ControlPoint>();

		/// <summary>
		/// List with the Elements of the <see cref="Model"/>.
		/// </summary>
		

		public IEnumerable<IElement> EnumerateElements() => ElementsDictionary.Values;
		public IEnumerable<INode> EnumerateNodes() => ControlPointsDictionary.Values;
		public IEnumerable<ISubdomain> EnumerateSubdomains() => PatchesDictionary.Values;

		public IElement GetElement(int elementID) => ElementsDictionary[elementID];
		public INode GetNode(int nodeID) => ControlPointsDictionary[nodeID];
		public ISubdomain GetSubdomain(int subdomainID) => PatchesDictionary[subdomainID];

		public int NumElements => ElementsDictionary.Count;
		public int NumNodes => ControlPointsDictionary.Count;
		public int NumSubdomains => PatchesDictionary.Count;

		public IDofSerializer DofSerializer { get; set; } = new StandardDofSerializer();

		/// <summary>
		/// List with the Elements of the <see cref="Model"/>.
		/// </summary>
		public IList<Element> Elements => ElementsDictionary.Values.ToList();

		/// <summary>
		/// Dictionary with the Elements of the <see cref="Model"/>.
		/// </summary>
		public Dictionary<int, Element> ElementsDictionary { get; } = new Dictionary<int, Element>();

		/// <summary>
		/// <see cref="IGlobalFreeDofOrdering"/> of the degrees of freedom of the <see cref="Model"/>.
		/// </summary>
		public IGlobalFreeDofOrdering GlobalDofOrdering { get; set; }
		//public IGlobalFreeDofOrdering GlobalDofOrdering
		//{
		//	get => _globalDofOrdering;
		//	set
		//	{
		//		_globalDofOrdering = value;
		//		foreach (var patch in Patches)
		//		{
		//			patch.FreeDofOrdering = GlobalDofOrdering.SubdomainDofOrderings[patch];
		//		}
		//	}
		//}

		private readonly List<PenaltyDofPair> penaltyBC = new List<PenaltyDofPair>();

		public void AddPenaltyConstrainedDofPair(PenaltyDofPair penaltyDofPair)
		{
			penaltyBC.Add(penaltyDofPair);
			var id = ElementsDictionary.Keys.Last() + 1;
			penaltyDofPair.ID = id;
			Element element = new PenaltyDofPair(penaltyDofPair.FirstPenaltyDof, penaltyDofPair.SecondPenaltyDof)
			{
				ID = id,
				Patch = PatchesDictionary[0],
				ElementType = penaltyDofPair,
			};

			PatchesDictionary[0].Elements.Add(element);
			ElementsDictionary.Add(id, element);
		}

		/// <summary>
		/// List containing the loads applied to the the <see cref="Model"/>.
		/// </summary>
		public IList<Load> Loads { get; private set; } = new List<Load>();

		/// <summary>
		/// List of <see cref="IMassAccelerationHistoryLoad"/> applied to the <see cref="Model"/>.
		/// </summary>
		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } =
			new List<IMassAccelerationHistoryLoad>();

		/// <summary>
		/// Return an <see cref="IReadOnlyList{ControlPoint}"/> with the Control Points of the <see cref="Model"/> as <see cref="INode"/>.
		/// </summary>
		//IReadOnlyList<INode> IModel.Nodes => ControlPointsDictionary.Values.ToList();

		/// <summary>
		/// Number of interfaces between patches.
		/// </summary>
		public int NumberOfInterfaces { get; set; }

		/// <summary>
		/// Number of patches of the <see cref="Model"/>.
		/// </summary>
		public int NumberOfPatches { get; set; }

		/// <summary>
		/// List with the patches of the model.
		/// </summary>
		public IList<Patch> Patches => PatchesDictionary.Values.ToList();

		/// <summary>
		/// Dictionary with the patches of the model.
		/// </summary>
		public Dictionary<int, Patch> PatchesDictionary { get; } = new Dictionary<int, Patch>();

		/// <summary>
		/// Dictionary with the patches of the model returned as <see cref="ISubdomain"/>.
		/// </summary>
		//IReadOnlyList<ISubdomain> IStructuralModel.Subdomains => PatchesDictionary.Values.ToList();

		/// <summary>
		/// List of time dependent loads added to the <see cref="Model"/>.
		/// </summary>
		public IList<ITimeDependentNodalLoad> TimeDependentNodalLoads { get; private set; } =
			new List<ITimeDependentNodalLoad>();

		/// <summary>
		/// Assigns nodal loads of the <see cref="Model"/>.
		/// </summary>
		/// <param name="distributeNodalLoads"><inheritdoc cref="NodalLoadsToSubdomainsDistributor"/></param>
		//public void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		//{
		//	foreach (var patch in PatchesDictionary.Values) patch.Forces.Clear();
		//	AssignNodalLoads(distributeNodalLoads);
		//	//Add possible penalty forces
		//}

		/// <summary>
		/// Assigns mass acceleration loads of the time step to the <see cref="Model"/>.
		/// </summary>
		/// <param name="timeStep">An <see cref="int"/> denoting the number of the time step.</param>
		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			throw new NotImplementedException();
		}

		/// <summary>
		/// Assigns nodal loads of the <see cref="Model"/>.
		/// </summary>
		/// <param name="distributeNodalLoads"><inheritdoc cref="NodalLoadsToSubdomainsDistributor"/></param>
		//public void AssignNodalLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		//{
  //          var globalNodalLoads = new Table<INode, IDofType, double>();
  //          foreach (var load in Loads)
  //          {
  //              var loadTable = load.CalculateLoad();
  //              foreach ((INode node, IDofType dof, double load) tuple in loadTable)
  //              {
  //                  if (tuple.node.Constraints.Any(x => x.DOF == tuple.dof)) continue;
  //                  if (globalNodalLoads.Contains(tuple.node, tuple.dof))
  //                  {
  //                      globalNodalLoads[tuple.node, tuple.dof] += tuple.load;
  //                  }
  //                  else
  //                  {
  //                      globalNodalLoads.TryAdd(tuple.node, tuple.dof, tuple.load);
  //                  }

  //              }
  //          }

  //          Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
  //          foreach (var idSubdomainLoads in subdomainNodalLoads)
  //          {
  //              PatchesDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
  //          }

		//	//var globalNodalLoads = new Table<INode, IDofType, double>();
		//	//foreach (var load in Loads) globalNodalLoads.TryAdd(load.Node, load.DOF, load.Amount);

		//	//var subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
		//	//foreach (var idSubdomainLoads in subdomainNodalLoads)
		//	//{
		//	//	PatchesDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
		//	//}
		//}

		/// <summary>
		/// Assigns mass acceleration loads of the time step to the <see cref="Model"/>.
		/// </summary>
		/// <param name="timeStep">An <see cref="int"/> denoting the number of the time step.</param>
		//public void AssignTimeDependentNodalLoads(int timeStep, NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		//{
		//	var globalNodalLoads = new Table<INode, IDofType, double>();
		//	foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
		//	{
		//		globalNodalLoads.TryAdd(load.Node, load.DOF, load.GetLoadAmount(timeStep));
		//	}

		//	Dictionary<int, SparseVector> subdomainNodalLoads = distributeNodalLoads(globalNodalLoads);
		//	foreach (var idSubdomainLoads in subdomainNodalLoads)
		//	{
		//		PatchesDictionary[idSubdomainLoads.Key].Forces.AddIntoThis(idSubdomainLoads.Value);
		//	}
		//}

		public void ApplyLoads()
		{
			foreach(var subdomain in PatchesDictionary.Values) subdomain.Forces.Clear();
			//AssignElementMassLoads();
			//AssignMassAccelerationLoads();
		}

		public void ApplyMassAccelerationHistoryLoads(int timeStep)
		{
			//if (MassAccelerationHistoryLoads.Count > 0)
			//{
			//	List<MassAccelerationLoad> m = new List<MassAccelerationLoad>(MassAccelerationHistoryLoads.Count);
			//	foreach (IMassAccelerationHistoryLoad l in MassAccelerationHistoryLoads)
			//	{
			//		m.Add(new MassAccelerationLoad() { Amount = l[timeStep], DOF = l.DOF });
			//	}

			//	foreach (Subdomain subdomain in SubdomainsDictionary.Values)
			//	{
			//		foreach (Element element in subdomain.Elements.Values)
			//		{
			//			double[] accelerationForces = element.ElementType.CalculateAccelerationForces(element, m);
			//			subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, accelerationForces, subdomain.Forces);
			//		}
			//	}
			//}

			//foreach (ElementMassAccelerationHistoryLoad load in ElementMassAccelerationHistoryLoads)
			//{
			//	MassAccelerationLoad hl = new MassAccelerationLoad()
			//	{
			//		Amount = load.HistoryLoad[timeStep] * 564000000,
			//		DOF = load.HistoryLoad.DOF
			//	};
			//	Element element = load.Element;
			//	ISubdomain subdomain = element.Subdomain;
			//	var accelerationForces = element.ElementType.CalculateAccelerationForces(
			//		load.Element, (new MassAccelerationLoad[] { hl }).ToList());
			//	GlobalDofOrdering.GetSubdomainDofOrdering(subdomain).AddVectorElementToSubdomain(element, accelerationForces,
			//		subdomain.Forces);
			//}
		}
		/// <summary>
		/// Clear the <see cref="Model"/>.
		/// </summary>
		public void Clear()
		{
			Loads.Clear();
			PatchesDictionary.Clear();
			ElementsDictionary.Clear();
			ControlPointsDictionary.Clear();
			_globalDofOrdering = null;
			Constraints.Clear();
			MassAccelerationHistoryLoads.Clear();
		}

		/// <summary>
		/// Interconnects Data Structures of the <see cref="Model"/>.
		/// </summary>
		public void ConnectDataStructures()
		{

			BuildInterconnectionData();
			AssignConstraints();
			AssignNodalLoadsToSubdomains();
			RemoveInactiveTransientNodalLoads();
		}

		private void AssignNodalLoadsToSubdomains()
		{
			// Remove inactive loads added by the user
			var activeLoads = new List<Load>(Loads.Count);
			foreach (Load load in Loads)
			{
				bool isConstrained = Constraints.Contains(load.Node, load.DOF);
				if (!isConstrained) activeLoads.Add(load);
			}
			Loads = activeLoads;

			// Assign the rest to their subdomains without scaling them. That will be done later by the analyzer and solver.
			foreach (Load load in Loads)
			{
				foreach (Patch subdomain in load.Node.SubdomainsDictionary.Values) subdomain.NodalLoads.Add(load);
			}
		}

		private void RemoveInactiveTransientNodalLoads()
		{
			var activeLoadsDynamic = new List<ITimeDependentNodalLoad>(TimeDependentNodalLoads.Count);
			foreach (ITimeDependentNodalLoad load in TimeDependentNodalLoads)
			{
				bool isConstrained = Constraints.Contains(load.Node, load.DOF);
				if (!isConstrained) activeLoadsDynamic.Add(load);
			}
			TimeDependentNodalLoads = activeLoadsDynamic;
		}

		private void AssignConstraints()
		{
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values)
			{
				if (controlPoint.Constraints == null) continue;
				foreach (Constraint constraint in controlPoint.Constraints)
					Constraints[controlPoint, constraint.DOF] = constraint.Amount;
			}

			foreach (Patch patch in PatchesDictionary.Values) patch.ExtractConstraintsFromGlobal(Constraints);
		}

		private void BuildElementDictionaryOfEachControlPoint()
		{
			foreach (var element in ElementsDictionary.Values)
			{
				foreach (var controlPoint in element.ControlPoints)
					controlPoint.ElementsDictionary.Add(element.ID, element);
			}
		}

		private void BuildInterconnectionData()
		{
			BuildPatchOfEachElement();
			BuildElementDictionaryOfEachControlPoint();
			foreach (var controlPoint in ControlPointsDictionary.Values) controlPoint.BuildPatchesDictionary();
		}

		private void BuildPatchOfEachElement()
		{
			foreach (var patch in PatchesDictionary.Values)
			{
				foreach (var element in patch.Elements)
				{
					element.Patch = patch;
					element.Model = this;
				}
			}
		}

	}
}
