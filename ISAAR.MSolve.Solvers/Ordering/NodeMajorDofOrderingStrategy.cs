﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: benchmark this against simple ordering + node major reordering
namespace ISAAR.MSolve.Solvers.Ordering
{
    /// <summary>
    /// Free dofs are assigned global / subdomain indices in a node major fashion: The dofs of the first node are 
    /// numbered, then the dofs of the second node, etc. Constrained dofs are ignored.
    /// </summary>
    public class NodeMajorDofOrderingStrategy : IFreeDofOrderingStrategy
    {
        public (int numGlobalFreeDofs, DofTable globalFreeDofs) OrderGlobalDofs(IStructuralModel model)
            => OrderFreeDofsOfElementSet(model.Elements, model.Nodes, model.Constraints);

        public (int numSubdomainFreeDofs, DofTable subdomainFreeDofs) OrderSubdomainDofs(ISubdomain subdomain)
            => OrderFreeDofsOfElementSet(subdomain.Elements, subdomain.Nodes, subdomain.Constraints);

        // Copied from the methods used by Subdomain and Model previously.
        private static (int numFreeDofs, DofTable freeDofs) OrderFreeDofsOfElementSet(IEnumerable<IElement> elements,
            IEnumerable<INode> sortedNodes, Table<INode, IDofType, double> constraints)
        {
            int totalDOFs = 0;
            Dictionary<int, List<IDofType>> nodalDOFTypesDictionary = new Dictionary<int, List<IDofType>>(); //TODO: use Set instead of List
            foreach (IElement element in elements)
            {
                var elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
                for (int i = 0; i < elementNodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(elementNodes[i].ID))
                        nodalDOFTypesDictionary.Add(elementNodes[i].ID, new List<IDofType>());
                    nodalDOFTypesDictionary[elementNodes[i].ID].AddRange(element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element)[i]);
                }
            }
            foreach (IElement element in elements)
            {
                for (int i = 0; i < element.Nodes.Count; i++)
                {
                    if (!nodalDOFTypesDictionary.ContainsKey(element.Nodes[i].ID))
                        nodalDOFTypesDictionary.Add(element.Nodes[i].ID, new List<IDofType>());
                    nodalDOFTypesDictionary[element.Nodes[i].ID].AddRange(element.ElementType.DofEnumerator.GetDofTypesForDofEnumeration(element)[i]);
                }
            }
            //foreach (IElement element in elements)
            //{
            //    var elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            //    for (int i = 0; i < elementNodes.Count; i++)
            //    {
            //        if (!nodalDOFTypesDictionary.ContainsKey(elementNodes[i].ID))
            //            nodalDOFTypesDictionary.Add(elementNodes[i].ID, new List<IDofType>());
            //        nodalDOFTypesDictionary[elementNodes[i].ID].AddRange(element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element)[i]);                    
            //        ///nodalDOFTypesDictionary[elementNodes[i].ID].AddRange(elementNodes[i].);
            //    }
            //}

            var freeDofs = new DofTable();
            foreach (INode node in sortedNodes)
            {
                //List<DOFType> dofTypes = new List<DOFType>();
                //foreach (Element element in node.ElementsDictionary.Values)
                //{
                //    if (elementsDictionary.ContainsKey(element.ID))
                //    {
                //        foreach (DOFType dof in element.ElementType.DOFTypes)
                //            dofTypes.Add(dof);
                //    }
                //}

                if (node.ID == 1129)
                {
                    string breakpoint = "@here";
                }

                Dictionary<IDofType, int> dofsDictionary = new Dictionary<IDofType, int>();
                //foreach (DOFType dofType in dofTypes.Distinct<DOFType>())
                foreach (IDofType dofType in nodalDOFTypesDictionary[node.ID].Distinct())
                {
                    int dofID = 0;
                    #region removeMaria
                    //foreach (DOFType constraint in node.Constraints)
                    //{
                    //    if (constraint == dofType)
                    //    {
                    //        dofID = -1;
                    //        break;
                    //    }
                    //}
                    #endregion

                    foreach (var constraint in node.Constraints) //TODO: access the constraints from the subdomain
                    {
                        if (constraint.DOF == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }

                    //var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                    ////if (node.EmbeddedInElement != null && node.EmbeddedInElement.ElementType.GetDOFTypes(null)
                    ////    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    ////    dofID = -1;
                    //if (embeddedNode != null && embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null)
                    //    .SelectMany(d => d).Count(d => d == dofType) > 0)
                    //    dofID = -1;

                    if (dofID == 0)
                    {
                        dofID = totalDOFs;
                        freeDofs[node, dofType] = dofID;
                        totalDOFs++;
                    }
                }
            }

            return (totalDOFs, freeDofs);
        }
    }
}
