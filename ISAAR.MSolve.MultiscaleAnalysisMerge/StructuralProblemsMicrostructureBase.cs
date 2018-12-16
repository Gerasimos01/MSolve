using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;

using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.MultiscaleAnalysis.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Providers;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.MultiscaleAnalysis
{
    public abstract class StructuralProblemsMicrostructureBase
    {
        public int SolverData { get; set; }

        public virtual Dictionary<int,Element> GetBoundaryFiniteElementsDictionary(Model model, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Element> boundaryElements = new Dictionary<int, Element>();

            foreach(Element element in model.Elements)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                for (int j = 0; j < elementDOFTypes.Count; j++)
                {
                    INode elementNode = matrixAssemblyNodes[j];
                    if (boundaryNodes.ContainsKey(elementNode.ID))
                    {
                        containsBoundaryNode = true;
                        break;
                    }
                }

                if (containsBoundaryNode)
                {
                    boundaryElements.Add(element.ID, element);
                }

            }

            return boundaryElements;

        }

        public virtual Dictionary<int, Element> GetBoundaryFiniteElementsDictionary(Subdomain subdomain, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Element> subdomainBoundaryElements = new Dictionary<int, Element>();

            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                bool containsBoundaryNode = false;

                var elementDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var matrixAssemblyNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                for (int j = 0; j < elementDOFTypes.Count; j++)
                {
                    INode elementNode = matrixAssemblyNodes[j];
                    if (boundaryNodes.ContainsKey(elementNode.ID))
                    {
                        containsBoundaryNode = true;
                        break;
                    }
                }

                if (containsBoundaryNode)
                {
                    subdomainBoundaryElements.Add(element.ID, element);
                }

            }

            return subdomainBoundaryElements;

        }

        public virtual Dictionary<int, Dictionary<int, Element>> GetSubdomainsBoundaryFiniteElementsDictionaries(Model model, Dictionary<int, Node> boundaryNodes)
        {
            Dictionary<int, Dictionary<int, Element>> subdomainsBoundaryElements = new Dictionary<int, Dictionary<int, Element>>();

            foreach (Subdomain subdomain in model.Subdomains)
            {
                Dictionary<int, Element> subdBoundaryElements = GetBoundaryFiniteElementsDictionary(subdomain, boundaryNodes);
                subdomainsBoundaryElements.Add(subdomain.ID, subdBoundaryElements);
            }

            return subdomainsBoundaryElements;
            
        }


        public virtual Dictionary<int, ILinearSystem> CreateNecessaryLinearSystems(Model model)
        {
            var linearSystems = new Dictionary<int, ILinearSystem>();
            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)" tou opoiu ta stoixeia ginontai access kai me ID
            {
                linearSystems.Add(subdomain.ID, new SkylineLinearSystem(subdomain.ID, subdomain.Forces));// prosoxh sto Id twn subdomain
            }

            return linearSystems;
        }

        public virtual ISolver GetAppropriateSolver(Dictionary<int, ILinearSystem> linearSystems)
        {
            //TODO: use solver data to create the chosen ISolver
            if (linearSystems.Keys.Count==1)
            {
                var solver = new SolverSkyline(linearSystems[1]); //linearSystems.ElementAt(0);
                return solver;
            }
            else
            {
                throw new NotImplementedException();
            }            
        }

        public virtual (NewtonRaphsonNonLinearAnalyzerDevelop, ProblemStructural,ElementStructuralStiffnessProvider, SubdomainGlobalMapping[]) AnalyzeMicrostructure(Model model, Dictionary<int, ILinearSystem> linearSystems, ISolver solver,
            int increments, int MaxIterations, int IterationsForMatrixRebuild, Dictionary<int, Dictionary<DOFType, double>> totalPrescribedBoundaryDisplacements,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Node> boundaryNodes, Dictionary<int, Vector> uInitialFreeDOFDisplacementsPerSubdomain)
        {
            #region Creation of nessesary analyzers for NRNLAnalyzer
            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            int totalSubdomains = model.Subdomains.Count;
            var linearSystemsArray = new ILinearSystem[totalSubdomains];
            var subdomainUpdaters = new NonLinearSubdomainUpdaterWithInitialConditions[totalSubdomains];
            var subdomainMappers = new SubdomainGlobalMapping[totalSubdomains]; int counter = 0;
            foreach (Subdomain subdomain in model.Subdomains)//TODO : or else "in model.SubdomainsDictionary.Values)"
            {
                linearSystemsArray[counter] = linearSystems[subdomain.ID];
                subdomainUpdaters[counter] = new NonLinearSubdomainUpdaterWithInitialConditions(subdomain);
                subdomainMappers[counter] = new SubdomainGlobalMapping(subdomain);
                counter++;
            }


            ElementStructuralStiffnessProvider elementProvider = new ElementStructuralStiffnessProvider();
            Dictionary<int, EquivalentContributionsAssebler> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            //equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler(model.SubdomainsDictionary[1], elementProvider));
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                equivalentContributionsAssemblers.Add(subdomain.ID, new EquivalentContributionsAssebler(subdomain, elementProvider));
            }
            #endregion

            #region Creation of Microstructure analyzer (NRNLdevelop temporarilly). 
            NewtonRaphsonNonLinearAnalyzerDevelop microAnalyzer = new NewtonRaphsonNonLinearAnalyzerDevelop(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers,
                provider, increments, model.TotalDOFs, uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalPrescribedBoundaryDisplacements, equivalentContributionsAssemblers);
            microAnalyzer.SetMaxIterations = MaxIterations;
            microAnalyzer.SetIterationsForMatrixRebuild = IterationsForMatrixRebuild;
            #endregion

            #region solution and update ------------->THA MPEI ENTOS KLASHS: of free converged displacements vectors;
            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, microAnalyzer, linearSystems);
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion

            return (microAnalyzer,provider,elementProvider,subdomainMappers);
        }
    }
}
