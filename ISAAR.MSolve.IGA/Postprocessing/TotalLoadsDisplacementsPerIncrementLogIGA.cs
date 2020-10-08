using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Utilities;

//TODO: This class should only extract data. How to output them (print in .txt, .xlsx, etc) should be done by different or
//      child classes
//TODO: Should this matrix write the results periodically (e.g. when a buffer fills up) instead of each increment?
//TODO: Extend it to extract data for many dofs simultaneously. It should also handles subdomains itself.
//TODO: If the analyzer ends abruptly (e.g. snap-through point in load control), that should be written here.
//TODO: Perhaps the file should not be opened and closed at each increment. Instead it should stay open, but then it should be 
//      disposed properly.
namespace ISAAR.MSolve.IGA.Postprocessing
{
    /// <summary>
    /// This does not work if the requested node belongs to an element that contains embedded elements.
    /// </summary>
    public class TotalLoadsDisplacementsPerIncrementLogIGA: ITotalLoadsDisplacementsPerIncrementLog
    {
        private readonly Patch subdomain;//.
        //private readonly Dictionary<INode, HashSet<DOFType>> monitorDofs = new Dictionary<INode, HashSet<DOFType>>();
        private readonly bool isMonitorDofFree;
        private readonly ControlPoint monitorNode;
        private readonly IDofType monitorDof;
        private readonly string outputFile;

        /// <summary>
        /// In case of displacement control, where there is a prescribed displacement at the monitored dof, we can only
        /// access the applied displacement which is scaled to 1/loadSteps at the beginning and remains constant during the iterations.
        /// Therefore we need to keep track of the previous displacements. This only happens for constrained dofs. 
        /// For free dofs this field is not used.
        /// </summary>
        private double currentTotalDisplacement = 0.0;
        //TODO: It should not be stored at all. Instead we should be able to access the total prescribed displacement from the analyzer
        
        public TotalLoadsDisplacementsPerIncrementLogIGA(Patch subdomain, int maxIncrementsExpected,
            ControlPoint monitorNode, IDofType monitorDof, string outputFile)
        {
            this.subdomain = subdomain;
            this.monitorNode = monitorNode;
            this.monitorDof = monitorDof;

            isMonitorDofFree = true;
            foreach (Constraint constraint in monitorNode.Constraints) //TODO: use LINQ instead of this
            {
                if (constraint.DOF == monitorDof)
                {
                    isMonitorDofFree = false;
                    break;
                }
            }

            this.outputFile = outputFile;
        }

        /// <summary>
        /// Writes the header.
        /// </summary>
        public void Initialize()
        {
            // If all subdomains use the same file, then we need to open it in append mode. 
            //TODO: Also that will not work in parallel for many subdomains.
            using (var writer = new StreamWriter(outputFile, false)) // do not append, since this is a new analysis
            {
                // Header
                writer.Write("Increment, Iteration, ResidualNorm");
                writer.Write($", Total displacement (Node {monitorNode.ID} - dof {monitorDof})");
                writer.WriteLine($", Total internal force (Node {monitorNode.ID} - dof {monitorDof})");
            }

            //if (!subdomain.Nodes.ContainsKey(monitorNode.ID)) throw new ArgumentException(
            //    "The requested monitor node does not belong to this subdomain");
        }

        /// <summary>
        /// This also writes to the output file.
        /// </summary>
        /// <param name="totalDisplacements">
        /// The total displacements (start till current iteration of current increment) of the subdomain.
        /// </param>
        /// <param name="totalInternalForces">
        /// The total internal right hand side forces (start till current iteration of current increment) of the subdomain.
        /// </param>
        public void LogTotalDataForIncrement(int incrementNumber, int currentIterationNumber, double errorNorm,
            IVectorView totalDisplacements, IVectorView totalInternalForces)
        {
            double displacement, force;
            bool isFreeDof = subdomain.FreeDofOrdering.FreeDofs.TryGetValue(monitorNode, monitorDof, out int subdomainDofIdx); //TODO: Should this be cached?
            if (isFreeDof)
            {
                displacement = totalDisplacements[subdomainDofIdx];
                force = totalInternalForces[subdomainDofIdx];
            }
            else
            {
                currentTotalDisplacement += subdomain.Constraints[monitorNode, monitorDof];
                displacement = currentTotalDisplacement;
                force = CalculateForceAt(monitorNode, monitorDof, totalDisplacements); //TODO: find out exactly what happens with the sign
            }

            // If all subdomains use the same file, then we need to open it in append mode. 
            //TODO: Also that will not work in parallel for many subdomains.
            using (var writer = new StreamWriter(outputFile, true)) // append mode to continue from previous increment
            {
                writer.Write($"{incrementNumber}, {currentIterationNumber}, {errorNorm}");
                writer.WriteLine($", {displacement}, {force}");
            }
        }

        internal double CalculateForceAt(ControlPoint node, IDofType dofType, IVectorView totalDisplacements)
        {
            double totalForce = 0.0;

            foreach (Element element in node.ElementsDictionary.Values)
            {
                // It is possible that one of the elements at this node does not engage this dof type, in which case -1 will be returned.
                // We will not have any contribution from them. If none of the elements engage this dof type, the total force will always be 0.
                int monitorDofIdx = FindLocalDofIndex(element, node, dofType);
                if (monitorDofIdx == -1) continue;

                //TODO: if an element has embedded elements, then we must also take into account their forces.
                double[] totalElementDisplacements = subdomain.CalculateElementDisplacements(element, totalDisplacements);
                double[] elementForces = element.ElementType.CalculateForcesForLogging(element, totalElementDisplacements);

                totalForce += elementForces[monitorDofIdx];
            }

            return totalForce;
        }

        /// <summary>
        /// Returns -1 if the element does not engage the requested <see cref="IDofType"/>
        /// </summary>
        private int FindLocalDofIndex(Element element, ControlPoint node, IDofType dofType)
        {
            int localNodeIdx = element.ControlNodes.IndexOf(node);
            Debug.Assert(localNodeIdx != -1, "The element does not contain this node.");
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            int localDofIdx = elementDofs[localNodeIdx].FindFirstIndex(dofType);
            int multNum = elementDofs[localNodeIdx].Count;
            int dofIdx = multNum * (localNodeIdx + 1) - (localDofIdx + 1);
            return dofIdx;
        }
    }
}
