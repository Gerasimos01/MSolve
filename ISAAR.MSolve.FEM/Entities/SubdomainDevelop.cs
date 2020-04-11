﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: remove code that calculates rhs vector components (nodal loads, constraints, etc). It should be moved to dedicated 
//      classes like EquivalentLoadAssembler, so that it can be reused between subdomains of different projects (FEM, IGA, XFEM).
//TODO: same for multiscale
namespace ISAAR.MSolve.FEM.Entities
{
    public class SubdomainDevelop : ISubdomain
    {
        public SubdomainDevelop(int id)
        {
            this.ID = id;
        }

        public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

        public Dictionary<int, Element> Elements { get; } = new Dictionary<int, Element>();

        //public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

        public int ID { get; }

        public List<Load> NodalLoads { get; } = new List<Load>();

        //TODO: I would prefer a Dictionary for fast random access, but a lot of tests break. Also having them sorted results in
        //      better orderings for structured meshes, although that can be implemented as a reordering or by ordering them 
        //      with LINQ before starting the ordering procedure
        public SortedDictionary<int, Node> Nodes { get; } = new SortedDictionary<int, Node>();

        public int NumElements => Elements.Count;
        public int NumNodalLoads => NodalLoads.Count;
        public int NumNodes => Nodes.Count;

        public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }
        public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

        public Vector Forces { get; set; } //TODO: this doesn't belong here

        //public bool MaterialsModified
        //{
        //    get
        //    {
        //        bool modified = false;
        //        foreach (Element element in elementsDictionary.Values)
        //            if (element.ElementType.MaterialModified)
        //            {
        //                modified = true;
        //                break;
        //            }
        //        return modified;
        //    }
        //}
        public bool StiffnessModified { get; set; } = true; // At first it is modified
        public bool ConnectivityModified { get; set; } = true; // At first it is modified

        //TODO: This belongs somewhere else (e.g. in EquivalentLoadsAssembler). It is not the Subdomain's job to calculate loading vectors
        public double[] CalculateElementDisplacements(Element element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses()
        {
            foreach (Element element in Elements.Values) element.ElementType.ClearMaterialStresses();
        }

        public void ConnectDataStructures()
        {
            DefineNodesFromElements();

            foreach (Element element in Elements.Values)
            {
                foreach (Node node in element.Nodes) node.ElementsDictionary[element.ID] = element;
            }
        }

        public void DefineNodesFromElements()
        {
            Nodes.Clear();
            foreach (Element element in Elements.Values)
            {
                foreach (Node node in element.Nodes) Nodes[node.ID] = node; // It may already be contained
            }

            //foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
            //    EmbeddedNodes.Add(e);
        }

        public IEnumerable<IElement> EnumerateElements() => Elements.Values;
        public IEnumerable<INodalLoad> EnumerateNodalLoads() => NodalLoads;
        public IEnumerable<INode> EnumerateNodes() => Nodes.Values;

        //TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
        //      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
        //      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
        //      It is too easy to access the wrong instance of the constraint. 
        //TODO: perhaps this code should be in Model. Assigning global data to each subdomain is usually done by the model.
        public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
        {
            //TODO: perhaps it is more efficient to traverse the global constraints instead of the subdomain's nodes, provided
            //      the latter are stored as a set. 
            //TODO: the next could be a Table method: Table.KeepDataOfRows(IEnumerable<TRow> rows)
            foreach (INode node in EnumerateNodes())
            {
                bool isNodeConstrained = globalConstraints.TryGetDataOfRow(node,
                    out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
                if (isNodeConstrained)
                {
                    foreach (var dofDisplacementPair in constraintsOfNode)
                    {
                        Constraints[node, dofDisplacementPair.Key] = dofDisplacementPair.Value;
                    }
                }
            }

            // This is probably faster but assumes that nodes store their prescribed displacements, which I hate.
            //foreach (Node node in Nodes)
            //{
            //    if (node.Constraints == null) continue;
            //    foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            //}
        }

        public IElement GetElement(int elementID) => Elements[elementID];
        public INode GetNode(int nodeID) => Nodes[nodeID];

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in Elements.Values)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

                //TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
                double[] localSolution = CalculateElementDisplacements(element, solution);
                double[] localdSolution = CalculateElementDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.StiffnessModified = true;
                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void CalculateStressesOnly(IVectorView solution, IVectorView dSolution)
        {
            foreach (Element element in Elements.Values)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

                //TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
                double[] localSolution = CalculateElementDisplacements(element, solution);
                double[] localdSolution = CalculateElementDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.StiffnessModified = true;
            }
            
        }

        public IVector CalculateRHSonly(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in Elements.Values)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

                //TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
                double[] localSolution = CalculateElementDisplacements(element, solution);
                double[] localdSolution = CalculateElementDisplacements(element, dSolution);

                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.StiffnessModified = false;
            foreach (Element element in Elements.Values) element.ElementType.ResetMaterialModified();
        }

        public void SaveMaterialState()
        {
            foreach (Element element in Elements.Values) element.ElementType.SaveMaterialState();
        }

        //TODO: I am against modifying the constraints table of the subdomain. Instead the analyzer should keep a constraint
        //      displacements vector at global/subdomain scale and modify that.
        public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

        //ADDED1
        public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, Node> boundaryNodes,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            ////prosthiki print

            //ekteleseis_counter2 += 1;
            //string counter_data = ekteleseis_counter2.ToString();
            //string path = string.Format(string2, counter_data);
            ////solution.WriteToFile(path);
            //double[] solution_data = new double[solution.Length];
            //solution_data = solution.CopyToArray();
            //WriteToFileVector(solution_data, path);

            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in Elements.Values)
            {
                var localSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
                ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                var localdSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
                ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(element, localdSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.StiffnessModified = true;
                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public double[] GetLocalVectorFromGlobalWithoutPrescribedDisplacements(Element element, IVectorView globalDisplacementVector)
        {
            double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            return elementNodalDisplacements;
        }
        

        public void ImposePrescribedDisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, Node> boundaryNodes,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {

            var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];
                int nodalDofsNumber = elementDOFTypes[j].Count;
                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                {
                    Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    int positionOfDofInNode = 0;
                    foreach (IDofType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDofInNode += 1;
                    }
                }
                iElementMatrixColumn += nodalDofsNumber;
            }

        }

        public void ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(Element element, double[] localSolution, Dictionary<int, Node> boundaryNodes,
            Dictionary<int, Dictionary<IDofType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<IDofType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {

            var elementDOFTypes = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];
                int nodalDofsNumber = elementDOFTypes[j].Count;
                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                {
                    Dictionary<IDofType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<IDofType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    int positionOfDofInNode = 0;
                    foreach (IDofType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            localSolution[iElementMatrixColumn + positionOfDofInNode] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                            // 1) den vazoume mono (1/increments) alla (nIncrement/increments) dioti metaxu aftwn twn nIncrements den exei mesolavhsei save sta material ths mikroklimakas
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDofInNode += 1;
                    }
                }
                iElementMatrixColumn += nodalDofsNumber;
            }

        }
    }
}